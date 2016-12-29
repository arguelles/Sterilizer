#!/bin/sh

function get_http_proxy(){
	if [ ! -z "$_get_http_proxy" ]; then
		return
	fi

	# Let's try to find a proxy...
	export OSG_SQUID_LOCATION=${OSG_SQUID_LOCATION:-UNAVAILABLE}
	if [ "$OSG_SQUID_LOCATION" != UNAVAILABLE ]; then
		export http_proxy=$OSG_SQUID_LOCATION
	fi

	# Try various proxies that may be reachable from CHTC
	if [ -z "$http_proxy" ]; then
		export http_proxy=osg-ss-se.chtc.wisc.edu:3128
		wget --spider -q http://convey.icecube.wisc.edu || unset http_proxy
	fi

	if [ -z "$http_proxy" ]; then
		export http_proxy=frontier01.hep.wisc.edu:3128
		wget --spider -q http://convey.icecube.wisc.edu || unset http_proxy
	fi

	if [ -z "$http_proxy" ]; then
		export http_proxy=frontier02.hep.wisc.edu:3128
		wget --spider -q http://convey.icecube.wisc.edu || unset http_proxy
	fi
	export _get_http_proxy=1
}

function fetch_globus(){
	if [ ! -z "$_fetch_globus" ]; then
		return
	fi
	# fetch globus toolkit if it's not already on the node
	if [ -d $CVMFS_ROOT/globus ]; then
		_GLOBUS_ROOT=$CVMFS_ROOT/globus
	else
		wget -q -U 'ICDiffuse(Wget)' --user icecube --password skua http://icecube:skua@x2100.icecube.wisc.edu/downloads/globus.tar.gz
		ERR=$?
		if [ $ERR -ne 0 -o ! -f globus.tar.gz ]; then
			echo "Failed to fetch globus tarball with wget" 1>&2
			exit 1
		fi
		tar xzf globus.tar.gz
		if [ ! -d globus ]; then
			echo "Failed to correctly unpack globus tarball" 1>&2
			exit 1
		fi
		_GLOBUS_ROOT=`pwd`/globus
	fi
	export X509_CERT_DIR=$_GLOBUS_ROOT/certificates
	export LD_LIBRARY_PATH=$_GLOBUS_ROOT/lib:$LD_LIBRARY_PATH
	export PATH=$_GLOBUS_ROOT/bin/:$PATH
	export _fetch_globus=1
}

DATASET=""
DATANAME=""
TARBALL="osg_diffuse.tar.xz"
DOWNLOAD_URL="gsiftp://gridftp.icecube.wisc.edu/data/user/cweaver/Diffuse/IC79/DiffuseAnalysis/LikelihoodFit/osg/"
PROG_DIR="osg_diffuse"
UPLOAD_URL="http://misc.hallsofchaos.net/Diffuse/upload.php"

CVMFS_ROOT=/cvmfs/oasis.opensciencegrid.org/icecube

# Handle arguments
while true; do
	if [ "$1" = "--dataset" ]; then
		shift
		DATASET=$1
		shift
	elif [ "$1" = "--name" ]; then
		shift
		DATANAME=$1
		shift
	elif [ "$1" = "--tarball" ]; then
		shift
		TARBALL=$1
		shift
	elif [ "$1" = "--input_url" ]; then
		shift
		DOWNLOAD_URL=$1
		shift
	elif [ "$1" = "--output_url" ]; then
		shift
		UPLOAD_URL=$1
		shift
	else
		break
	fi
done

if [ -z "$DATASET" ]; then
	echo "Dataset name must be specified with --dataset" 1>&2
	exit 1
fi
if [ -z "$DATANAME" ]; then
	echo "Data file name must be specified with --name" 1>&2
	exit 1
fi

if [ -n "$OSG_GRID" ]; then
	source $OSG_GRID/setup.sh
fi

get_http_proxy
fetch_globus

# Fetch tarball
if [ ! -f $TARBALL ]; then
	globus-url-copy -nodcau ${DOWNLOAD_URL}${TARBALL} file://`pwd`/$TARBALL
	ERR=$?

	if [ $ERR -ne 0 -o ! -f $TARBALL ]; then
		echo "host: "`hostname`" proxy: "$http_proxy" globus: "$_GLOBUS_ROOT 1>&2
		echo "Problem downloading tarball" 1>&2
		exit 1
	fi
fi

# Unpack tarball
tar xJf $TARBALL
ERR=$?

if [ $ERR -ne 0 -o ! -d ${PROG_DIR} ]; then
	echo "Problem unpacking tarball" 1>&2
	exit 2
fi

# Prepare environment
cd $PROG_DIR
ERR=$?

if [ $ERR -ne 0 ]; then
	echo "Failed to enter program directory" 1>&2
	exit 3
fi

PROG_DIR=`pwd`
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${PROG_DIR}/lib"
export NEWNUFLUX_DATADIR="${PROG_DIR}/flux_paramterizations/"

# Run!
./IC79_diffuse "$@" --readCompact -o=output.txt
ERR=$?

if [ $ERR -ne 0 ]; then
	echo "Calculation failed" 1>&2
	exit 4
fi

# Report results
if `echo "$UPLOAD_URL" | grep -q '^http://'`; then
	# send results to upload.php using wget
	echo -n "store=true&dataset=" > results.txt
	echo -n "$DATASET" | ./urlencode >> results.txt
	echo -n "&name=" >> results.txt
	echo -n "$DATANAME" | ./urlencode >> results.txt
	echo -n "&data=" >> results.txt
	cat output.txt | ./urlencode >> results.txt
	wget -q -U 'ICDiffuse(Wget)' --user=icecube --password=skua --post-file=results.txt $UPLOAD_URL
	ERR=$?

	if [ $ERR -ne 0 ]; then
		echo "HTTP upload failed" 1>&2
		exit 5
	fi
elif `echo "$UPLOAD_URL" | grep -q '^gsiftp://'`; then
	# upload with GridFTP
	cd ..
	globus-url-copy file://${PROG_DIR}/output.txt ${UPLOAD_URL}/${DATASET}/${DATANAME}
	ERR=$?

	if [ $ERR -ne 0 ]; then
		echo "GridFTP upload failed" 1>&2
		exit 5
	fi
elif `echo "$UPLOAD_URL" | grep -q '^file://'`; then
	# copy back with cp
	UPLOAD_URL=`echo "$UPLOAD_URL" | sed 's|^file://||'`
	cp output.txt ${UPLOAD_URL}/${DATASET}/${DATANAME}

	if [ $ERR -ne 0 ]; then
		echo "cp upload failed" 1>&2
		exit 5
	fi
else
	echo "Unrecognized upload protocol" 1>&2
	exit 6
fi
