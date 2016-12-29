#ifndef LF_TABLEIO_H
#define LF_TABLEIO_H

#include <cassert>
#include <string>
#include <map>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/mpl/string.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

#include <hdf5.h>
#include <hdf5_hl.h>

template <int N>
constexpr char at(char const(&s)[N], int i){
	return(i >= N ? '\0' : s[i]);
}

template <class S, char C, bool EOS>
struct push_back_if_c :
boost::mpl::push_back<S, boost::mpl::char_<C>>
{};

template <class S, char C>
struct push_back_if_c<S, C, true> : S {};

#ifndef CTS_STRING_MAX_LENGTH
#define CTS_STRING_MAX_LENGTH 48
#endif

#define CTS_PRE(z, n, u) push_back_if_c<
#define CTS_POST(z, n, s) , at(s, n), (n >= sizeof(s))>::type

#define CTS(s)                                   \
BOOST_PP_REPEAT(CTS_STRING_MAX_LENGTH, CTS_PRE, ~)    \
boost::mpl::string<>                       \
BOOST_PP_REPEAT(CTS_STRING_MAX_LENGTH, CTS_POST, s)

//Type is the type of the field, Name is an MPL sequence from which the
//field name can be derived by application of mpl::c_str
template<typename Type, typename Name>
struct field{
	typedef Type type;
	typedef Name name; //note that this is a type, not a value!
	static std::string getName(){
		return(boost::mpl::c_str<Name>::value);
	}
};

namespace detail{
	template<class... FieldTypes> struct hoboTuple;
	
	template<>
	struct hoboTuple<>{
		static const int nFields=0;
		
		static void getFieldSizes(size_t sizes[], unsigned int counter){}
		static void getFieldOffsets(size_t offsets[], unsigned int counter, size_t base){}
		static void getFieldNames(std::string names[], unsigned int counter){}
	};
	
	template<unsigned int index, class FieldType, class... OtherFields>
	struct FieldType_{
		typedef typename FieldType_<index-1,OtherFields...>::type type;
	};
	
	template<class FieldType, class... OtherFields>
	struct FieldType_<0, FieldType, OtherFields...>{
		typedef typename FieldType::type type;
	};
	
	template<unsigned int index, class FieldType, class... OtherFields>
	struct getImpl{
		static typename FieldType_<index,FieldType,OtherFields...>::type get(const hoboTuple<FieldType,OtherFields...>& tup){
			return(getImpl<index-1,OtherFields...>::get(tup.suffix));
		}
		static typename FieldType_<index,FieldType,OtherFields...>::type& get(hoboTuple<FieldType,OtherFields...>& tup){
			return(getImpl<index-1,OtherFields...>::get(tup.suffix));
		}
	};
	
	template<class FieldType, class... OtherFields>
	struct getImpl<0,FieldType,OtherFields...>{
		static typename FieldType_<0,FieldType,OtherFields...>::type get(const hoboTuple<FieldType,OtherFields...>& tup){
			return(tup.field);
		}
		static typename FieldType_<0,FieldType,OtherFields...>::type& get(hoboTuple<FieldType,OtherFields...>& tup){
			return(tup.field);
		}
	};
	
	template<class FieldType, class... OtherFields>
	struct hoboTuple<FieldType,OtherFields...>{
		typedef hoboTuple<OtherFields...> SuffixType;
		
		typename FieldType::type field;
		SuffixType suffix;
		
		static const int nFields=SuffixType::nFields+1;
		
		hoboTuple() = default;
		
		template<typename Arg, typename... ArgTypes>
		hoboTuple(Arg arg, ArgTypes... args):field(arg),suffix(args...){}
		
		static void getFieldSizes(size_t sizes[], unsigned int counter){
			sizes[counter++]=sizeof(typename FieldType::type);
			SuffixType::getFieldSizes(sizes,counter);
		}
		
		static void getFieldOffsets(size_t offsets[], unsigned int counter, size_t base){
			typedef hoboTuple<FieldType,OtherFields...> ThisType;
			offsets[counter++]=base+offsetof(ThisType,field);
			SuffixType::getFieldOffsets(offsets,counter,base+offsetof(ThisType,suffix));
		}
		
		static void getFieldNames(std::string names[], unsigned int counter){
			names[counter++]=FieldType::getName();
			SuffixType::getFieldNames(names,counter);
		}
	};
	
	template<unsigned int index_, typename type_>
	struct IndexForNameHolder{
		static constexpr unsigned int index=index_;
		typedef type_ type;
	};
	
	//fwd decl
	template<unsigned int index_, typename Name, typename FieldType, typename... OtherFieldTypes>
	struct IndexForNameImpl;
	
	//base case
	template<unsigned int index_, typename Name, typename FieldType>
	struct IndexForNameImpl<index_,Name,FieldType>{
	private:
		static_assert(std::is_same<Name,typename FieldType::name>::value, "Name does not name a field in the current tuple");
		using underlying_type = typename std::conditional<std::is_same<Name,typename FieldType::name>::value,
		IndexForNameHolder<index_,typename FieldType::type>,
		void
		>::type;
	public:
		static constexpr unsigned int index=underlying_type::index;
		using type = typename underlying_type::type;
	};
	
	//recursive case
	template<unsigned int index_, typename Name, typename FieldType, typename... OtherFieldTypes>
	struct IndexForNameImpl{
	private:
		static_assert(std::is_same<Name,typename FieldType::name>::value || sizeof...(OtherFieldTypes)>0, "Name does not name a field in the current tuple");
		using underlying_type = typename std::conditional<std::is_same<Name,typename FieldType::name>::value,
		IndexForNameHolder<index_,typename FieldType::type>,
		IndexForNameImpl<index_+1,Name,OtherFieldTypes...>
		>::type;
	public:
		static constexpr unsigned int index=underlying_type::index;
		using type = typename underlying_type::type;
	};
	
	template<typename Name, typename... FieldTypes>
	struct IndexForName{
	private:
		using underlying_type = IndexForNameImpl<0,Name,FieldTypes...>;
	public:
		static constexpr unsigned int index=underlying_type::index;
		using type = typename underlying_type::type;
	};
	
}

struct RecordID{
	uint32_t run;
	uint32_t event;
	uint32_t subEvent;
    int subEventStream;
	bool exists;
	
	bool operator==(const RecordID& other) const{
		return(run==other.run && event==other.event && subEventStream==other.subEventStream && subEvent==other.subEvent);
	}
	
	bool operator<(const RecordID& other) const{
		if(run<other.run)
			return(true);
		if(run>other.run)
			return(false);
		if(event<other.event)
			return(true);
		if(event>other.event)
			return(false);
		if(subEventStream<other.subEventStream)
			return(true);
		if(subEventStream>other.subEventStream)
			return(false);
		if(subEvent<other.subEvent)
			return(true);
		return(false);
	}
};

template<class... FieldTypes>
struct TableRow{
	typedef detail::hoboTuple<FieldTypes...> FieldsType;
	
	RecordID id;
	FieldsType fields;
	
	static const int baseFields=5;
	static const int nFields=baseFields+FieldsType::nFields;
	
	TableRow() = default;
	
	template<typename... ArgTypes>
	TableRow(ArgTypes... args):fields(args...){}
	
	static void getFieldSizes(size_t sizes[]){
		sizes[0]=sizeof(uint32_t);
		sizes[1]=sizeof(uint32_t);
		sizes[2]=sizeof(uint32_t);
		sizes[3]=sizeof(int);
		sizes[4]=sizeof(bool);
		FieldsType::getFieldSizes(sizes,baseFields);
	}
	
	static void getFieldOffsets(size_t offsets[]){
		typedef TableRow<FieldTypes...> ThisType;
		offsets[0]=offsetof(ThisType,id.run);
		offsets[1]=offsetof(ThisType,id.event);
		offsets[2]=offsetof(ThisType,id.subEvent);
		offsets[3]=offsetof(ThisType,id.subEventStream);
		offsets[4]=offsetof(ThisType,id.exists);
		FieldsType::getFieldOffsets(offsets,baseFields,offsetof(ThisType,fields));
	}
	
	static void getFieldNames(std::string names[nFields]){
		names[0]="Run";
		names[1]="Event";
		names[2]="SubEvent";
		names[3]="SubEventStream";
		names[4]="exists";
		FieldsType::getFieldNames(names,baseFields);
	}
	
	template<unsigned int index>
	typename detail::FieldType_<index,FieldTypes...>::type get() const{
		static_assert(index<=FieldsType::nFields,"Attempt to access out of range column");
		return(detail::getImpl<index,FieldTypes...>::get(fields));
	}
	
	template<unsigned int index>
	typename detail::FieldType_<index,FieldTypes...>::type& get(){
		static_assert(index<=FieldsType::nFields,"Attempt to access out of range column");
		return(detail::getImpl<index,FieldTypes...>::get(fields));
	}
	
	template<typename Name, typename Info=detail::IndexForName<Name,FieldTypes...>>
	const typename Info::type get() const{
		static_assert(Info::index<=FieldsType::nFields,"Attempt to access out of range column");
		return(detail::getImpl<Info::index,FieldTypes...>::get(fields));
	}
	
	template<typename Name, typename Info=detail::IndexForName<Name,FieldTypes...>>
	typename Info::type get(){
		static_assert(Info::index<=FieldsType::nFields,"Attempt to access out of range column");
		return(detail::getImpl<Info::index,FieldTypes...>::get(fields));
	}
};

struct H5File{
	hid_t file;
	H5File(std::string path, unsigned int flags = H5F_ACC_RDONLY, hid_t fapl_id = H5P_DEFAULT):
	file(H5Fopen(path.c_str(), flags, fapl_id)){}
	~H5File(){
		H5Fclose(file);
	}
	operator hid_t() const{ //implicit conversion
		return(file);
	}
	operator bool() const{ //conversion to bool, tests if file is valid
		return(file>=0);
	}
};

namespace{
	//this is supposed to live in readTableBlock, but ATM using a local struct for a deleter seems to make clang crash
	struct nameDeleter{
		unsigned int size;
		nameDeleter(unsigned int s):size(s){}
		void operator()(char** obj){
			for(unsigned int i=0; i<size; i++)
				delete[] (obj[i]);
			delete[] obj;
		}
	};
}

//read up to 1 MB of data at a time
template<typename T, typename Event, typename Callback>
void readTable(hid_t tableLoc, const std::string& tableName, std::map<RecordID,Event>& data, Callback callback){
	hsize_t nFields,nRecords;
	H5TBget_table_info(tableLoc,tableName.c_str(),&nFields,&nRecords);
	//TODO: check for errors
	assert(nFields>=T::nFields);
	std::unique_ptr<char*[],nameDeleter> availableFieldNames(new char*[nFields],nameDeleter(nFields));
	for(unsigned int i=0; i<nFields; i++)
		availableFieldNames[i]=new char[255];
	std::unique_ptr<size_t[]> availableFieldSizes(new size_t[nFields]);
	std::unique_ptr<size_t[]> availableFieldOffsets(new size_t[nFields]);
	size_t tableRowSize;
	H5TBget_field_info(tableLoc,tableName.c_str(),availableFieldNames.get(),availableFieldSizes.get(),availableFieldOffsets.get(),&tableRowSize);
	//TODO: check for errors
	
	std::unique_ptr<std::string[]> requestedFieldNames(new std::string[T::nFields]);
	T::getFieldNames(requestedFieldNames.get());
	
	int fieldIndices[T::nFields];
	size_t fieldSizes[T::nFields], fieldOffsets[T::nFields];
	T::getFieldSizes(fieldSizes);
	T::getFieldOffsets(fieldOffsets);
	unsigned int lastIndex=0;
	
	for(unsigned int i=0; i<T::nFields; i++){
		for(; lastIndex<nFields; lastIndex++){
			if(availableFieldNames[lastIndex]==requestedFieldNames[i])
				break;
		}
		if(lastIndex==nFields)
			throw std::runtime_error("Field '"+requestedFieldNames[i]+
									 "' either does not exist in table '"+tableName+"' or was requested out of order");
		fieldIndices[i]=lastIndex;
		if(availableFieldSizes[lastIndex]!=fieldSizes[i])
			throw std::runtime_error("Field '"+requestedFieldNames[i]+
									 "' size does not match: target size is "+boost::lexical_cast<std::string>(fieldSizes[i])+
									 " but size in table is "+boost::lexical_cast<std::string>(availableFieldSizes[lastIndex]));
	}
	
	const unsigned int maxBlockSize=1<<20; //1 MB
	unsigned int recordsPerBlock=maxBlockSize/sizeof(T);
	if(!recordsPerBlock)
		recordsPerBlock=1;
	std::vector<T> intermediate(recordsPerBlock);
	
	hsize_t recordsRead=0;
	while(recordsRead<nRecords){
		hsize_t toRead=recordsPerBlock;
		if(recordsRead+toRead>nRecords)
			toRead=nRecords-recordsRead;
		intermediate.resize(toRead);
		T* dataPtr=&intermediate[0];
		if(H5TBread_fields_index(tableLoc, tableName.c_str(), T::nFields, fieldIndices, recordsRead, toRead, sizeof(T),  fieldOffsets, fieldSizes, dataPtr)<0)
			throw std::runtime_error("Read error");
		for(T& sourceItem : intermediate){
			Event& targetItem=data[sourceItem.id];
			callback(sourceItem,targetItem);
		}
		recordsRead+=toRead;
	}
}

#endif