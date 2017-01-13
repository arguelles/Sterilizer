#ifndef timer
#define timer


struct timer{
private:
  std::chrono::high_resolution_clock::time_point t1, t2;
  bool running=false;
  std::string name;
  friend std::ostream& operator<<(std::ostream& os, timer& t);
public:
timer():running(false){}
  explicit timer(std::string name):running(false),name(name){}

  ~timer(){
    if(running)
      std::cout << *this << std::endl;
  }

  std::string getName() const{
    return(name);
  }
  void setName(std::string name){
    this->name=name;
  }

  void start(){
    if(!running){
      t1 = std::chrono::high_resolution_clock::now();
      running=true;
    }
  }
  void start(std::string name){
    if(!running){
      this->name=name;
      t1 = std::chrono::high_resolution_clock::now();
      running=true;
    }
  }


  void stop(){
    if(running){
      t2 = std::chrono::high_resolution_clock::now();
      running=false;
    }
  }
  double get(){
    if(running)
      stop();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    return(time_span.count());
  }
};

std::ostream& operator<<(std::ostream& os, timer& t){
  if(quiet)
    return(os);
  if(!t.getName().empty())
    os << t.getName() << ": ";
  return(os << t.get() << " seconds");
}




#endif
