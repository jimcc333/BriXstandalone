#include "generaltools.h"

// returns true if data is increasing (cumulative sum)
bool DecreaseChecker(vector<float> &data) {

    const int length = data.size();

    if(length < 2) {return false;}

    for(int i = 1; i < length; i++) {
        if(data[i-1] > data[i]) {return true;}
    }

    return false;
}

void CumulativeAdd(vector<float> &data) {
    const int length = data.size();

    if(length < 2) {return;}

    for(int i = 1; i < length; i++) {
        data[i] += data[i-1];
    }
}

// Linear interpolation function
float Interpolate(const float y0, const float y1, const float x0, const float x1, const float x) {
    return y0 + (y1 - y0)*(x - x0)/(x1 - x0);
}

void UpdateStorage(std::vector<int> &decay_times){
    for(int i = 0; i < decay_times.size(); i++){decay_times[i]++;}
    return;
}
