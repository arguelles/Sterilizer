#ifndef LF_SPLINETABLE_H
#define LF_SPLINETABLE_H

#include <unistd.h>
#include <numeric>
#include <iostream>
#include <cmath>
#include <cassert>

#include <photospline/splinetable.h>
#include <photospline/bspline.h>
#include <photospline/glam.h>

#include <boost/lexical_cast.hpp>
#include <PhysTools/plottable_histograms.h>

//#include "plottingExtras.h"

//wrapper for sucky C code
class Splinetable{
public:
	splinetable data;
	int searchCenterSuccess;
	bool readFromFile;
	
	void computeStrides(){
		data.strides = malloc<unsigned long>(data.ndim);
		long arraysize;
		data.strides[data.ndim - 1] = arraysize = 1;
		for (int i = data.ndim-1; i >= 0; i--) {
			arraysize *= data.naxes[i];
			if (i > 0)
				data.strides[i-1] = arraysize;
		}
	}
	
	///type-safe wrapper for malloc
	template<typename T>
	static T* malloc(size_t count){
		return(static_cast<T*>(::malloc(count*sizeof(T))));
	}
	
	void clear(){
		//attempt to compensate for shenanigans
		if(!readFromFile){
			for(unsigned int i=0; i<data.ndim; i++){
				if(data.knots[i])
					data.knots[i]+=data.order[i];
			}
		}
		splinetable_free(&data);
	}
	
public:
	Splinetable():readFromFile(false){
		memset(&data, 0, sizeof(data));
	}
	Splinetable(const Splinetable& other){
		data.ndim=other.data.ndim;
		data.order=malloc<int>(data.ndim);
		std::copy_n(other.data.order,data.ndim,data.order);
		
		data.nknots=malloc<long>(data.ndim);
		std::copy_n(other.data.nknots,data.ndim,data.nknots);
		data.knots=malloc<double*>(data.ndim);
		for(unsigned int i=0; i<data.ndim; i++){
			data.knots[i]=malloc<double>(data.nknots[i]);
			std::copy_n(other.data.knots[i],data.nknots[i],data.knots[i]);
		}
		
		//TODO: copy extents
		data.extents=nullptr;
		//TODO: copy periods
		data.periods=nullptr;
		
		data.naxes=malloc<long>(data.ndim);
		std::copy_n(other.data.naxes,data.ndim,data.naxes);
		
		unsigned long nCoeffs=std::accumulate(data.naxes, data.naxes+data.ndim, 1UL, std::multiplies<unsigned long>());
		data.coefficients=malloc<float>(nCoeffs);
		std::copy_n(other.data.coefficients,nCoeffs,data.coefficients);
		
		computeStrides();
	}
	//construct by reading from a file
	Splinetable(std::string tablefile):readFromFile(true){
		readsplinefitstable(tablefile.c_str(), &data);
		assert(data.ndim>0);
	}

  void read(std::string tablefile){
    readFromFile=true;
    readsplinefitstable(tablefile.c_str(), &data);
    assert(data.ndim>0);
  }

	//disallow assignment
	Splinetable& operator=(const Splinetable&)=delete;
	
	~Splinetable(){
		clear();
	}
	
	size_t getDimension() const{ return(data.ndim); }
	
	//evaluate the spline at a point whose coordinates are in the array pointed to by x,
	//which is assumed to have data.ndim entries
	double operator()(const double* x){
		int centerBuffer[data.ndim];
		if((searchCenterSuccess=tablesearchcenters(&data, x, centerBuffer)) == 0)
			return(ndsplineeval(&data, x, centerBuffer, 0));
		return(0.0);
	}
	//same as evaluation operator, but compute the derivative at x with respect to dimension dim
	double derivative(const double* x, unsigned int dim){
		int centerBuffer[data.ndim];
		if((searchCenterSuccess=tablesearchcenters(&data, x, centerBuffer)) == 0)
			return(ndsplineeval(&data, x, centerBuffer, 1<<dim));
		return(0.0);
	}
	
	//write to a file, overwirting if the file already exists
	bool write(std::string tablefile){
		unlink(tablefile.c_str());
		return(writesplinefitstable(tablefile.c_str(), &data)==0);
	}
	
	template<typename HistType, typename TransformType=std::function<double(unsigned int,double)>, typename BinType=typename HistType::dataType>
	void fitHistogram(const HistType& h, std::vector<std::vector<double>> knots,
				 std::vector<unsigned int> splineOrder, std::vector<double> smoothing,
				 TransformType coordinateTransform=[](unsigned int dim, double d){ return(d); }){
		//need to have the right number of knot vectors
		assert(h.getDimensions()==knots.size());
		//need to have either one or a matching number of penatly strengths
		assert(smoothing.size()==1 || smoothing.size()==h.getDimensions());
		
		//throw out any existing data
		clear();
		
		//set dimensions
		data.ndim=h.getDimensions();
		//copy spline orders
		data.order=malloc<int>(data.ndim);
		std::copy(splineOrder.begin(),splineOrder.end(),&data.order[0]);
		//set penalty orders (currently fixed to 2)
		int* penaltyOrder=new int[data.ndim];
		std::fill_n(&penaltyOrder[0],data.ndim,2);
		//copy knot locations
		data.nknots=malloc<long>(data.ndim);
		data.knots=malloc<double*>(data.ndim);
		for(unsigned int i=0; i<knots.size(); i++){
			data.nknots[i]=knots[i].size();
			data.knots[i]=malloc<double>(data.nknots[i]);
			std::copy(knots[i].begin(), knots[i].end(), data.knots[i]);
		}
		
		auto binCheck=[](const BinType& bin)->bool{
			using EVT=phys_tools::ErrorValueTraits<BinType>;
			return(!std::isinf(EVT::value(bin)) && !std::isnan(EVT::value(bin))
				   && EVT::error(bin)!=0  && std::isnormal(EVT::error(bin)));
		};
		
		//first pass: figure out how many valid points we have
		ndsparse sparseData;
		sparseData.rows=0;
		sparseData.ndim=h.getDimensions();
		for(const auto bin : h){
			if(binCheck(bin))
				sparseData.rows++;
		}
		//allocate space
		sparseData.i=new int*[sparseData.ndim];
		sparseData.ranges=new int[sparseData.ndim];
		double** binCoordinates=new double*[sparseData.ndim];
		for(unsigned int i=0; i<h.getDimensions(); i++){
			sparseData.ranges[i]=h.getBinCount(i);
			sparseData.i[i]=new int[sparseData.rows];
			binCoordinates[i]=new double[sparseData.ranges[i]];
			for(unsigned int j=0; j<sparseData.ranges[i]; j++)
				binCoordinates[i][j]=coordinateTransform(i,h.getBinCenter(i,j));
		}
		//second pass: copy data to the ndsparse and compute weights
		sparseData.x=new double[sparseData.rows];
		double* weights=new double[sparseData.rows];
		{
			unsigned long i=0;
			for(auto it=h.cbegin(), end=h.cend(); it!=end; it++){
				if(binCheck(*it)){
					sparseData.x[i]=(double)(typename HistType::dataType)*it;
					for(int j=0; j<sparseData.ndim; j++)
						sparseData.i[j][i]=it.getCoordinate(j);
					weights[i]=1/((typename HistType::dataType)*it).error();
					assert(std::isnormal(weights[i]));
					i++; //must be done only for valid entries
				}
			}
		}
		
		cholmod_common cholmod_state;
		cholmod_l_start(&cholmod_state);
		//build penalty matrix
		cholmod_sparse* penalty;
		{
			std::unique_ptr<long[]> nsplines(new long[data.ndim]);
			unsigned long sidelen = 1;
			for(unsigned int i = 0; i < data.ndim; i++) {
				nsplines[i] = data.nknots[i] - data.order[i] - 1;
				sidelen *= nsplines[i];
			}
			penalty=cholmod_l_spzeros(sidelen, sidelen, 1, CHOLMOD_REAL, &cholmod_state);
			for(unsigned int i = 0; i < data.ndim; i++){
				penalty = add_penalty_term(nsplines.get(), data.knots[i], data.ndim,
										   i, data.order[i], penaltyOrder[i],
										   (smoothing.size()>1?smoothing[i]:smoothing[0]),
										   false/*no monodim*/, penalty, &cholmod_state);
			}
		}
		
		//do the fit
		glamfit_complex(&sparseData, weights, binCoordinates,
		                &data, data.order, penalty,
		                -1/*no monotonic dimension*/, 1/*not verbose*/,
		                &cholmod_state);
		cholmod_l_free_sparse(&penalty, &cholmod_state);
		cholmod_l_finish(&cholmod_state);
		
		//clean up buffers
		delete[] penaltyOrder;
		delete[] weights;
		delete[] sparseData.ranges;
		for(unsigned int i=0; i<sparseData.ndim; i++){
			delete[] binCoordinates[i];
			delete[] sparseData.i[i];
		}
		delete[] binCoordinates;
		delete[] sparseData.i;
		
		//fix up things which will be out of date
		computeStrides();
	}
	
	///make this table into a stacked interpolation of the given otehr tables
	///\param tables the tables to stack
	///\param coordinates the coordiantes in the stacking dimension of the tables to be stacked
	///\param stackOrder the order of the spline in the stacking dimension
	void stackTables(const std::vector<Splinetable*>& tables, std::vector<double> coordinates, int stackOrder=2){
		assert(!tables.empty());
		assert(tables.size()==coordinates.size());
		int inputDim=tables.front()->data.ndim;
		for(auto table : tables){
			assert(table->data.ndim == inputDim);
			for(unsigned int i=0; i<inputDim; i++)
				assert(table->data.order[i] && tables.front()->data.order[i]);
		}
		
		//throw out any existing data
		clear();
		
		//set dimensions
		data.ndim=inputDim+1;
		//copy/set spline orders
		data.order=malloc<int>(data.ndim);
		std::copy_n(tables.front()->data.order,inputDim,data.order);
		data.order[inputDim]=stackOrder;
		
		//copy/set knots
		data.nknots=malloc<long>(data.ndim);
		std::copy_n(tables.front()->data.nknots,inputDim,data.nknots);
		data.nknots[inputDim]=tables.size()+stackOrder+1;
		data.knots=malloc<double*>(data.ndim);
		//copy existing knots
		for(unsigned int i=0; i<inputDim; i++){
			data.knots[i]=malloc<double>(data.nknots[i]);
			std::copy_n(tables.front()->data.knots[i],data.nknots[i],data.knots[i]);
		}
		//figure out knots for the new dimension
		{
			data.knots[inputDim]=malloc<double>(data.nknots[inputDim]);
			double* lastKnots=data.knots[inputDim];
			//copy input positions
			std::copy_n(coordinates.begin(),tables.size(),lastKnots+stackOrder);
				
			//shift knots
			double knotShift=(stackOrder-1)*(lastKnots[stackOrder+tables.size()-1]-lastKnots[stackOrder])/(2*tables.size());
			for(unsigned int i=0; i<tables.size(); i++)
				lastKnots[stackOrder+i]+=knotShift;
			
			//add stackOrder padding knots before
			double knotStep=lastKnots[stackOrder+1]-lastKnots[stackOrder];
			for(int i=0; i<stackOrder; i++)
				lastKnots[i]=lastKnots[stackOrder]+(i-stackOrder)*knotStep;
			//add one padding knot after
			lastKnots[data.nknots[inputDim]-1]=2*lastKnots[data.nknots[inputDim]-2]-lastKnots[data.nknots[inputDim]-3];
		}
		
		//TODO: handle any extents or periods
		
		//set naxes
		data.naxes=malloc<long>(data.ndim);
		std::copy_n(tables.front()->data.naxes,inputDim,data.naxes);
		data.naxes[inputDim]=tables.size();
		
		//copy coefficients
		unsigned long nCoeffs=std::accumulate(data.naxes, data.naxes+data.ndim, 1UL, std::multiplies<unsigned long>());
		unsigned long nInputCoeffs=std::accumulate(data.naxes, data.naxes+data.ndim-1, 1UL, std::multiplies<unsigned long>());
		data.coefficients=malloc<float>(nCoeffs);
		unsigned int step=data.naxes[data.ndim-1];
		for(unsigned int i=0; i<tables.size(); i++){
			for(unsigned int j=0; j<nInputCoeffs; j++)
				data.coefficients[i+j*step]=tables[i]->data.coefficients[j];
		}
		
		computeStrides();
	}
};

#endif
