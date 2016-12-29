#include <functional>
#include <map>
#include <string>
#include <vector>
#include <iostream>

#include <boost/lexical_cast.hpp>

class OptionParser{
private:
	std::map<char,std::function<void(std::string)>> shortOptions;
	std::map<char,std::function<void()>> shortOptionsNoStore;
	std::map<std::string,std::function<void(std::string)>> longOptions;
	std::map<std::string,std::function<void()>> longOptionsNoStore;
	bool printedUsage;
	std::string usageMessage;
	
	void checkIdentifier(std::string ident){
		if(ident.empty())
			throw std::logic_error("Invalid option name: '': options may not be empty");
		if(ident.find('=')!=std::string::npos)
			throw std::logic_error("Invalid option name: '"+ident+"': options may not contain '='");
		if(ident.find('-')==0)
			throw std::logic_error("Invalid option name: '"+ident+"': options may not begin with '-'");
		//valid identifier, do nothing
	}
	
	bool optionKnown(char ident){
		return(shortOptions.count(ident) || shortOptionsNoStore.count(ident));
	}
	bool optionKnown(std::string ident){
		return(longOptions.count(ident) || longOptionsNoStore.count(ident));
	}
	
	template<typename IDType, typename DestType>
	void addOption(IDType ident, DestType& destination, std::map<IDType,std::function<void(std::string)>>& storeMap){
		if(optionKnown(ident))
			throw std::logic_error("Attempt to redefine option '"+boost::lexical_cast<std::string>(ident)+"'");
		storeMap.emplace(ident,
						 [ident,&destination](std::string optData)->void{
							 std::istringstream ss(optData);
							 ss >> destination;
							 if(ss.fail()){
								 throw std::runtime_error("Failed to parse \""+optData+"\" as argument to '"
														  +boost::lexical_cast<std::string>(ident)+"' option");
							 }
						 });
	}
	template<typename IDType>
	void addOption(IDType ident, std::function<void()> action, std::map<IDType,std::function<void()>>& storeMap){
		if(optionKnown(ident))
			throw std::logic_error("Attempt to redefine option '"+boost::lexical_cast<std::string>(ident)+"'");
		storeMap.emplace(ident,action);
	}
	
	bool handleNextArg(const std::string& arg){
		static const auto npos=std::string::npos;
		if(arg.size()<2) //not an option, skip it
			return(false);
		if(arg[0]!='-') //not an option, skip it
			return(false);
		size_t startIdx=arg.find_first_not_of('-');
		if(startIdx>2) //not an option, skip it
			return(false);
		size_t endIdx=arg.find('=',startIdx);
		std::string opt=arg.substr(startIdx,(endIdx==npos?npos:endIdx-startIdx));
		
		if(opt.empty())
			throw std::runtime_error("Invalid option: '"+arg+"'");
		
		if((opt.size()==1 && startIdx>1) || (opt.size()>1 && startIdx==1))
			throw std::runtime_error("Malformed option: '"+arg+"' (wrong number of leading dashes)");
		
		std::string value;
		if(endIdx!=npos && endIdx!=arg.size()-1)
			value=arg.substr(endIdx+1);
		
		if(opt.size()==1){
			char optC=opt[0];
			if(shortOptions.count(optC)){
				if(endIdx==npos)
					throw std::runtime_error("Malformed option: '"+arg+"' (missing value after '=')");
				shortOptions.find(optC)->second(value);
			}
			else if(shortOptionsNoStore.count(optC)){
				if(endIdx!=npos)
					throw std::runtime_error("Malformed option: '"+arg+"' (no value expected for this flag)");
				shortOptionsNoStore.find(optC)->second();
			}
			else
				throw std::runtime_error("Unknown option: '"+arg+"'");
		}
		else{
			if(longOptions.count(opt)){
				if(endIdx==npos)
					throw std::runtime_error("Malformed option: '"+arg+"' (missing value after '=')");
				longOptions.find(opt)->second(value);
			}
			else if(longOptionsNoStore.count(opt)){
				if(endIdx!=npos)
					throw std::runtime_error("Malformed option: '"+arg+"' (no value expected for this flag)");
				longOptionsNoStore.find(opt)->second();
			}
			else
				throw std::runtime_error("Unknown option: '"+arg+"'");
		}
		return(true);
	}
	
	/*template<typename Iterator>
	bool handleNextArg(Iterator it, Iterator end){
		static const auto npos=std::string::npos;
		
		const std::string& arg=*it;
		
		if(arg.size()<2) //not an option, skip it
			return(false);
		if(arg[0]!='-') //not an option, skip it
			return(false);
		size_t startIdx=arg.find_first_not_of('-');
		if(startIdx>2) //not an option, skip it
			return(false);
		
		size_t endIdx=arg.find('=',startIdx);
		std::string opt=arg.substr(startIdx,(endIdx==npos?npos:endIdx-startIdx));
		
		if(opt.empty())
			throw std::runtime_error("Invalid option: '"+arg+"'");
		
		std::string value;
		if(endIdx!=npos && endIdx!=arg.size()-1)
			value=arg.substr(endIdx+1);
	}*/
	
	static std::string synonymList(const std::initializer_list<std::string>& list){
		std::ostringstream ss;
		for(auto begin=list.begin(), it=begin, end=list.end(); it!=end; it++){
			if(it!=begin)
				ss << ", ";
			if(it->size()==1)
				ss << '-';
			else
				ss << "--";
			ss << *it;
		}
		return(ss.str());
	}
	
public:
	explicit OptionParser(bool automaticHelp=true):printedUsage(false){
		if(automaticHelp)
			addOption({"h","?","help","usage"},
					  [this](){
						  std::cout << getUsage() << std::endl;
						  printedUsage=true;
					  },
					  "Print usage information.");
	}
	
	bool didPrintUsage() const{
		return(printedUsage);
	}
	
	template<typename T>
	void addOption(char ident, T& destination, std::string description){
		checkIdentifier(std::string(1,ident));
		addOption(ident,destination,shortOptions);
		std::ostringstream ss;
		ss  << " -" << ident << ": " << description << '\n';
		usageMessage+=ss.str();
	}
	void addOption(char ident, std::function<void()> action, std::string description){
		checkIdentifier(std::string(1,ident));
		addOption(ident,std::move(action),shortOptionsNoStore);
		std::ostringstream ss;
		ss  << " -" << ident << ": " << description << '\n';
		usageMessage+=ss.str();
	}
	template<typename T>
	void addOption(std::string ident, T& destination, std::string description){
		checkIdentifier(ident);
		addOption(ident,destination,longOptions);
		std::ostringstream ss;
		ss  << " --" << ident << ": " << description << '\n';
		usageMessage+=ss.str();
	}
	void addOption(std::string ident, std::function<void()> action, std::string description){
		checkIdentifier(ident);
		addOption(ident,std::move(action),longOptionsNoStore);
		std::ostringstream ss;
		ss  << " --" << ident << ": " << description << '\n';
		usageMessage+=ss.str();
	}
	template<typename T>
	void addOption(std::initializer_list<std::string> idents, T& destination, std::string description){
		for(auto ident : idents)
			checkIdentifier(ident);
		for(auto ident : idents){
			if(ident.size()==1)
				addOption(ident[0],destination,shortOptions);
			else
				addOption(ident,destination,longOptions);
		}
		std::ostringstream ss;
		ss << ' ' << synonymList(idents) << ": " << description << '\n';
		usageMessage+=ss.str();
	}
	void addOption(std::initializer_list<std::string> idents, std::function<void()> action, std::string description){
		for(auto ident : idents)
			checkIdentifier(ident);
		for(auto ident : idents){
			if(ident.size()==1)
				addOption(ident[0],action,shortOptionsNoStore);
			else
				addOption(ident,action,longOptionsNoStore);
		}
		std::ostringstream ss;
		ss << ' ' << synonymList(idents) << ": " << description << '\n';
		usageMessage+=ss.str();
	}
	
	template<typename Iterator>
	std::vector<std::string> parseArgs(Iterator argBegin, Iterator argEnd){
		std::vector<std::string> positionals;
		while(argBegin!=argEnd){
			std::string arg=*argBegin;
			if(!handleNextArg(arg))
				positionals.push_back(arg);
			argBegin++;
		}
		return(positionals);
	}
	std::vector<std::string> parseArgs(int argc, char* argv[]){
		return(parseArgs(argv,argv+argc));
	}
	
	std::string getUsage(){
		return(usageMessage);
	}
};