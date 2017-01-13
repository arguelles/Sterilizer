#include "SterileSearch.h"

void SterileSearch::LoadData(std::string dataPath){
		try{
			sample=loadExperimentalData(dataPath,UseBurnsample);
		} catch(std::exception& ex){
			std::cerr << "Problem loading experimental data: " << ex.what() << std::endl;
			return(1);
		}
		if(!quiet)
			std::cout << "Loaded " << sample.size() << " experimental events" << std::endl;
}

void SterileSearch::LoadCompactData(std::string filepath) {

}

marray<double,3> SterileSearch::GetDataDistribution(){
    marray<double,3> array {static_cast<size_t>(data_hist.getBinCount(2)),
                            static_cast<size_t>(data_hist.getBinCount(1)),
                            static_cast<size_t>(data_hist.getBinCount(0))};

    for(size_t iy=0; iy<data_hist.getBinCount(2); iy++){ // year
      for(size_t ic=0; ic<data_hist.getBinCount(1); ic++){ // zenith
        for(size_t ie=0; ie<data_hist.getBinCount(0); ie++){ // energy
          auto itc = static_cast<likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(data_hist(ie,ic,iy));
          array[iy][ic][ie] = itc.size();
        }
      }
    }
    return array;
}

marray<double,3> SterileSearch::GetExpectation(SterileNeutrinoParameters snp, std::vector<double> nuisance){
    MakeSimulationHistogram(snp,nuisance);
    marray<double,3> array {static_cast<size_t>(sim_hist.getBinCount(2)),
                            static_cast<size_t>(sim_hist.getBinCount(1)),
                            static_cast<size_t>(sim_hist.getBinCount(0))};

    auto weighter = DFWM(nuisance);
    for(size_t iy=0; iy<sim_hist.getBinCount(2); iy++){ // year
      for(size_t ic=0; ic<sim_hist.getBinCount(1); ic++){ // zenith
        for(size_t ie=0; ie<sim_hist.getBinCount(0); ie++){ // energy
          auto itc = static_cast<likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(sim_hist(ie,ic,iy));
          double expectation=0;
          for(auto event : itc.entries()){
            expectation+=weighter(event);
          }
          array[iy][ic][ie] = expectation;
        }
      }
    }
    return array;
}

marray<double,3> SterileSearch::GetRealization(SterileNeutrinoParameters snp, std::vector<double> nuisance){

}
