#ifndef PLOTTINGEXTRAS_H
#define PLOTTINGEXTRAS_H

#include <iostream>
#include <cmath>
#include <cassert>

#include <boost/lexical_cast.hpp>
#include <PhysTools/plottable_histograms.h>

template<typename T>
struct interval{
	T min,max;
	
	bool contains(const T& x) const{
		return(x>=min && x<=max);
	}
};

interval<unsigned int> fcPoissonAcceptanceInterval(double mean, double confidence, double offset=0.0);
interval<double> fcPoissonConfidenceInterval(unsigned int obs, double confidence, double offset=0.0, double tol=1e-4);

struct fcErrorValue{
private:
	double value, errMin, errMax;
	
	void recomputeErrors(){
		if(value<100){
			interval<double> errRange=fcPoissonConfidenceInterval(value, .68, 0.0, 1e-3);
			errMin=errRange.min;
			errMax=errRange.max;
		}
		else{
			double r=sqrt(value);
			errMin=value-r;
			errMax=value+r;
		}
	}
public:
	fcErrorValue():value(0){}
	
	fcErrorValue(double d):value(d){
		recomputeErrors();
	}
	
	fcErrorValue& operator+=(const double& d){
		value+=d;
		recomputeErrors();
		return(*this);
	}
	fcErrorValue& operator+=(const fcErrorValue& other){
		value+=other.value;
		recomputeErrors();
		return(*this);
	}
	fcErrorValue& operator*=(double scale){
		value*=scale;
		recomputeErrors();
		return(*this);
	}
	fcErrorValue operator+(double a) const{
		return(fcErrorValue(value+a));
	}
	fcErrorValue operator+(const fcErrorValue& other) const{
		return(fcErrorValue(value+other.value));
	}
	fcErrorValue operator*(double scale) const{
		return(fcErrorValue(value*scale));
	}
	
	operator double() const{ return(value); }
	double errorMin() const{ return(errMin); }
	double errorMax() const{ return(errMax); }	
};

fcErrorValue operator*(double scale, const fcErrorValue& v);
std::ostream& operator<<(std::ostream& os, const fcErrorValue& v);

namespace phys_tools{
	namespace histograms{
		namespace detail{
			template<>
			class histogramTraits<fcErrorValue>{
			public:
				using amount=amount_type<double>;
				constexpr static bool enable_automatic_amount_handling=true;
				static void defaultData(fcErrorValue* data, unsigned int count){/* Nothing to do */}
				static double unit(){ return 1.0; }
			};
		}
	}
}
	
struct sqErrorValue{
private:
	double value, value2;
public:
	sqErrorValue():value(0),value2(0){}
	
	sqErrorValue(double d):value(d),value2(d*d){}
	
	sqErrorValue(double d, double e):value(d),value2(e*e){}
	
	void setError(double e){
		value2=(e*e);
	}
	
	sqErrorValue& operator+=(const double& d){
		value+=d;
		value2+=d*d;
		return(*this);
	}
	sqErrorValue& operator+=(const sqErrorValue& other){
		value+=other.value;
		value2+=other.value2;
		return(*this);
	}
	sqErrorValue& operator*=(double scale){
		value*=scale;
		value2*=scale*scale;
		return(*this);
	}
	sqErrorValue& operator/=(double scale){
		value/=scale;
		value2/=scale*scale;
		return(*this);
	}
	//danger! this assumes that the two values have uncorrelated errors
	sqErrorValue& operator*=(const sqErrorValue& other){
		double p=value*other.value;
		double rn=value2/(value*value);
		double rd=other.value2/(other.value*other.value);
		value=p;
		value2=p*p*(rn+rd);
		return(*this);
	}
	//danger! this assumes that the two values have uncorrelated errors
	sqErrorValue& operator/=(const sqErrorValue& other){
		double r1=value/other.value;
		double rn=value2/(value*value);
		double rd=other.value2/(other.value*other.value);
		value=r1;
		value2=r1*r1*(rn+rd);
		return(*this);
	}
	sqErrorValue operator+(double a) const{
		return(sqErrorValue(*this)+=a);
	}
	sqErrorValue operator+(const sqErrorValue& other) const{
		return(sqErrorValue(*this)+=other);
	}
	sqErrorValue operator*(double scale) const{
		return(sqErrorValue(*this)*=scale);
	}
	sqErrorValue operator/(double scale) const{
		return(sqErrorValue(*this)/=scale);
	}
	//danger! this assumes that the two values have uncorrelated errors
	sqErrorValue operator*(const sqErrorValue& other) const{
		return(sqErrorValue(*this)*=other);
	}
	//danger! this assumes that the two values have uncorrelated errors
	sqErrorValue operator/(const sqErrorValue& other) const{
		return(sqErrorValue(*this)/=other);
	}
	
	operator double() const{ return(value); }
	double error() const{ return(sqrt(value2)); }
	double errorMin() const{ return(value-sqrt(value2)); }
	double errorMax() const{ return(value+sqrt(value2)); }
};

sqErrorValue operator*(double scale, const sqErrorValue& v);
std::ostream& operator<<(std::ostream& os, const sqErrorValue& v);

namespace phys_tools{
	namespace histograms{
		namespace detail{
			template<>
			class histogramTraits<sqErrorValue>{
			public:
				using amount=amount_type<sqErrorValue>;
				constexpr static bool enable_automatic_amount_handling=true;
				static void defaultData(sqErrorValue* data, unsigned int count){/* Nothing to do */}
				static double unit(){ return 1.0; }
			};
		}
	}
	
	template<>
	class ErrorValueTraits<sqErrorValue>{
	public:
		static double value(const sqErrorValue& t){ return((double)t); }
		static double error(const sqErrorValue& t){ return(t.error()); }
	};
}

/*class PlottableFilledCurves : public phys_tools::gnuplot::Plottable{
public:
	//min and max _must_ have identical binning!
	phys_tools::histograms::histogram<1> min;
	phys_tools::histograms::histogram<1> max;
	
	PlottableFilledCurves(const std::string title, const phys_tools::histograms::histogram<1>& min, const phys_tools::histograms::histogram<1>& max):
	phys_tools::gnuplot::Plottable(title,phys_tools::gnuplot::Format::FILLED_CURVES),min(min),max(max){}
	
	PlottableFilledCurves(const std::string title, phys_tools::histograms::histogram<1>&& min, phys_tools::histograms::histogram<1>&& max):
	phys_tools::gnuplot::Plottable(title,phys_tools::gnuplot::Format::FILLED_CURVES),min(std::move(min)),max(std::move(max)){}
	
	virtual std::string get_data_code() const{
		unsigned int count=min.getBinCount(0);
		std::ostringstream data;
		for(unsigned int i=0; i<count; i++){
			double x=min.getBinEdge(0,i);
			data << x << '\t' << min(i) << '\t' << max(i) << '\n';
			x=min.getBinEdge(0,i)+min.getBinWidth(0,i);
			data << x << '\t' << min(i) << '\t' << max(i) << '\n';
		}
		data << 'e' << std::endl;
		return(data.str());
	}
	
	virtual std::string get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const{
		return(get_data_code());
	}
	
	virtual std::string get_plot_code () const{
        std::ostringstream plot_code;
        plot_code << "'-' title \"" << title_ << "\" ";
        plot_code << format_;
        return (plot_code.str());
    }
};*/
	
template<typename BinType>
class FilledErrorHist : public phys_tools::histograms::histogram<1,BinType>, public phys_tools::gnuplot::Plottable{
public:
	class Format : public phys_tools::gnuplot::Format{
	private:
		phys_tools::Value<phys_tools::gnuplot::ColorSpec>  fill_color_;
		phys_tools::gnuplot::Format::Style error_style_;
		bool showFill_, showErrors_;
	public:
		Format():phys_tools::gnuplot::Format(FILLED_CURVES),error_style_(LEFT_STEPS),showFill_(true),showErrors_(true){}
		Format(phys_tools::gnuplot::Format::Style s):phys_tools::gnuplot::Format(s),error_style_(LEFT_STEPS),showFill_(true),showErrors_(true){}
		Format(phys_tools::gnuplot::Format f):phys_tools::gnuplot::Format(f),error_style_(LEFT_STEPS),showFill_(true),showErrors_(true){}
		
		phys_tools::gnuplot::ColorSpec fill_color() const{ return(fill_color_); }
		Format& fill_color(phys_tools::gnuplot::ColorSpec fill_color){
			fill_color_=fill_color;
			return(*this);
		}
		
		phys_tools::gnuplot::Format::Style error_style() const{ return(error_style_); }
		Format& error_style(phys_tools::gnuplot::Format::Style es){
			error_style_=es;
			return(*this);
		}
		
		bool showFill() const{ return(showFill_); }
		Format& showFill(bool show){
			showFill_=show;
			return(*this);
		}
		
		bool showErrors() const{ return(showErrors_); }
		Format& showErrors(bool show){
			showErrors_=show;
			return(*this);
		}
		
		Format makeFillFormat() const{
			Format f(*this);
			if(style()==FILLED_CURVES && fill_color_.good())
				f.line_color(fill_color_);
			return(f);
		}
		Format makeErrorFormat() const{
			Format f(*this);
			f.style(error_style_).point_size(0);
			return(f);
		}
	};
private:
	Format format_;
public:
	using iterator=typename phys_tools::histograms::histogram<1,BinType>::iterator;
	using const_iterator=typename phys_tools::histograms::histogram<1,BinType>::const_iterator;

	template<typename AxisType, typename = typename std::is_convertible<typename std::add_lvalue_reference<AxisType>::type,phys_tools::histograms::axis&>::type>
	FilledErrorHist(const std::string& title,Format f,AxisType axis):
	Plottable(title),phys_tools::histograms::histogram<1,BinType>(axis),format_(f){}
	
	FilledErrorHist(phys_tools::histograms::histogram<1,BinType>&& other):phys_tools::histograms::histogram<1,BinType>(std::move(other)){}
	
	FilledErrorHist(const std::string& title, phys_tools::histograms::histogram<1,BinType>&& other):Plottable(title),phys_tools::histograms::histogram<1,BinType>(std::move(other)){}
	
	Format& format(){ return(format_); }
	const Format& format() const{ return(format_); }
	
	phys_tools::histograms::histogram<1> getCenter(){
		using namespace phys_tools::histograms;
		histogram<1> center;
		center.setAxis(0,this->getAxis(0)->copy());
		bool saved=this->getUseContentScaling();
		this->setUseContentScaling(false);
		center.setUseContentScaling(false);
		for(const_iterator it=this->cbegin(), end=this->cend(); it!=end; it++)
			center.add(it.getBinCenter(0),amount((double)*it));
		this->setUseContentScaling(saved);
		center.setUseContentScaling(saved);
		return(center);
	}
	
	phys_tools::histograms::histogram<1> getMin(){
		using namespace phys_tools::histograms;
		histogram<1> min;
		min.setAxis(0,this->getAxis(0)->copy());
		bool saved=this->getUseContentScaling();
		this->setUseContentScaling(false);
		min.setUseContentScaling(false);
		for(const_iterator it=this->cbegin(), end=this->cend(); it!=end; it++)
			min.add(it.getBinCenter(0),amount((*it).errorMin()));
		this->setUseContentScaling(saved);
		min.setUseContentScaling(saved);
		return(min);
	}
	
	phys_tools::histograms::histogram<1> getMax(){
		using namespace phys_tools::histograms;
		histogram<1> max;
		max.setAxis(0,this->getAxis(0)->copy());
		bool saved=this->getUseContentScaling();
		this->setUseContentScaling(false);
		max.setUseContentScaling(false);
		for(const_iterator it=this->cbegin(), end=this->cend(); it!=end; it++)
			max.add(it.getBinCenter(0),amount((*it).errorMax()));
		this->setUseContentScaling(saved);
		max.setUseContentScaling(saved);
		return(max);
	}
	
	virtual std::string get_plot_code() const{
        std::ostringstream plot_code;
		if(format_.style()==phys_tools::gnuplot::Format::FILLED_CURVES){
			if(format_.showFill())
				plot_code << "'-' title \"" << title_ << "\" " << format_.makeFillFormat();
			if(format_.showFill() && format_.showErrors())
				plot_code << ", ";
			if(format_.showErrors()){
				plot_code << "'-' notitle " << format_.makeErrorFormat();
				plot_code << ", '-' notitle " << format_.makeErrorFormat();
			}
		}
		else{
			if(format_.showFill())
				plot_code << "'-' title \"" << title_ << "\" " << format_.makeFillFormat();
			if(format_.showFill() && format_.showErrors())
				plot_code << ", ";
			if(format_.showErrors())
				plot_code << "'-' notitle " << format_.makeErrorFormat();
		}
        return (plot_code.str());
    }
	
private:
	std::string get_data_code_core(const double x_min=-std::numeric_limits<double>::infinity(),
									const double x_max=std::numeric_limits<double>::infinity(),
									const double y_min=-std::numeric_limits<double>::infinity(),
									const double y_max=std::numeric_limits<double>::infinity()) const{
		unsigned int count=this->getBinCount(0);
		std::ostringstream data;
		auto clamp=[y_min,y_max](double y){
			if(y<y_min)
				return(y_min);
			if(y>y_max)
				return(y_max);
			return(y);
		};
		
		if(format_.style()==phys_tools::gnuplot::Format::FILLED_CURVES){
			auto outputErrors=[&,this](int bin, int item){
				const double min=clamp((*this)(bin).errorMin());
				const double max=clamp((*this)(bin).errorMax());
				
				switch(item){
					case 0: data << min << '\t' << max; break;
					case 1: data << min; break;
					case 2: data << max; break;
					
				}
				data << '\n';
			};
			
			const bool doublePoints=format_.error_style()==phys_tools::gnuplot::Format::LEFT_STEPS;
			auto outputPoint=
			[&,this](int bin, int item){
				data << this->getBinEdge(0,bin) << '\t';
				outputErrors(bin,item);
				if(doublePoints){
					data << this->getBinEdge(0,bin)+this->getBinWidth(0,bin) << '\t';
					outputErrors(bin,item);
				}
			};
			if(format_.showFill()){
				for(unsigned int i=0; i<count; i++)
					outputPoint(i,0);
				data << 'e' << std::endl;
			}
			if(format_.showErrors()){
				for(int item=1; item<3; item++){
					for(unsigned int i=0; i<count; i++)
						outputPoint(i,item);
					data << 'e' << std::endl;
				}
			}
		}
		else{
			const auto cStyle=format_.style();
			auto outputData=[&](int item){
				for(unsigned int i=0; i<count; i++){
					double x=0.0;
					switch(cStyle){
						case phys_tools::gnuplot::Format::LEFT_STEPS:
							x=this->getBinEdge(0,i);
							break;
						case phys_tools::gnuplot::Format::RIGHT_STEPS:
							x=this->getBinEdge(0,i)+this->getBinWidth(0,i);
							break;
						case phys_tools::gnuplot::Format::CENTER_STEPS:
						case phys_tools::gnuplot::Format::LINES:
						case phys_tools::gnuplot::Format::POINTS:
						case phys_tools::gnuplot::Format::LINES_POINTS:
						case phys_tools::gnuplot::Format::IMPULSES:
						case phys_tools::gnuplot::Format::DOTS:
						case phys_tools::gnuplot::Format::ERROR_BARS:
						case phys_tools::gnuplot::Format::Y_ERROR_BARS:
							//case phys_tools::gnuplot::Format::XY_ERROR_BARS:
							x=this->getBinCenter(0,i);
							break;
						default:
							throw std::runtime_error("Unsupported drawing format: "+boost::lexical_cast<std::string>(format().style()));
					}
					data << x;
					if(item==0)
						data << '\t' << (double)(*this)(i);
					if(item)
						data << '\t' << clamp((double)(*this)(i)) << '\t' << clamp((*this)(i).errorMin()) << '\t' << clamp((*this)(i).errorMax());
					data << '\n';
				}
				data << 'e' << std::endl;
			};
			if(format_.showFill())
				outputData(0);
			if(format_.showErrors())
				outputData(1);
		}
		return(data.str());
	}
public:
	virtual std::string get_data_code() const{
		return(get_data_code_core());
	}
	
	virtual std::string get_data_code (const double& x_min,
                                       const double& x_max,
                                       const double& y_min,
                                       const double& y_max) const{
		return(get_data_code_core(x_min,x_max,y_min,y_max));
	}
};

#endif