#include "iostream"
#include "fstream"
#include "string"
#include <map>

using namespace std;
multimap< int, pair<int, int> > JsonMaps[2];    // <runNumber, pair<startLumiSec, endLumiSec> >
void checkChars(char *nameInput, string &nameOutput, int &status);

bool isGoodEvt(int runNumber,int LumiSec){
    bool isGoodEvt_ = false;
    multimap< int, pair<int, int> >::iterator JsonMapItr_;
    multimap< int, pair<int, int> > JsonMap(JsonMaps[1]);    // <runNumber, pair<startLumiSec, endLumiSec> > for tempotary
    if (runNumber<=163869) JsonMap = JsonMaps[0];

    JsonMapItr_ = JsonMap.find(runNumber);

    if(JsonMapItr_ == JsonMap.end()){
        isGoodEvt_ = false;
    }else if(JsonMapItr_ != JsonMap.end()) {
        while (JsonMapItr_!= JsonMap.end()){
            if(JsonMapItr_->first!=runNumber) break;
            if(LumiSec>=JsonMapItr_->second.first&&LumiSec<=JsonMapItr_->second.second){
                isGoodEvt_ = true;
            }
            if(isGoodEvt_) break;
            JsonMapItr_++;
        }
    }
    return isGoodEvt_;
}

void MakeJsonMap(string jsonfile){
    string jsonfiles[2] = {jsonfile,jsonfile};

    for(int jsonindx_=0;jsonindx_<2;jsonindx_++){
       ifstream JSON(jsonfiles[jsonindx_].c_str());

       if(!JSON) {
           std::cout<<"[ERROR] Can not found JSON file, "<<jsonfiles[jsonindx_]<<". Please check if the JSON file exists."<<std::endl;
           exit(0);
       }   


       char name[128];
       int runNumber = 0;
       int startLumiSec = 0;
       int endLumiSec = 0;
       while(!JSON.eof()){
           JSON >> name;

           string nameOutput;
           int status;
           checkChars(name,nameOutput,status);

           if(status==1){
               runNumber = atoi(nameOutput.c_str());
           }else if(status==2){
               startLumiSec = atoi(nameOutput.c_str());
           }else if(status==3){
               endLumiSec = atoi(nameOutput.c_str());
               JsonMaps[jsonindx_].insert( pair< int, pair<int, int> >(runNumber, pair<int, int>(startLumiSec,endLumiSec)));
           }
       }
       JSON.close();
    }
}

void checkChars(char *nameInput, string &nameOutput, int &status){
        for(unsigned int size_=0;size_<strlen(nameInput);size_++){
            if(nameInput[0]=='{'||nameInput[0]=='"'){
                status = 1; // for run number
            }else if(nameInput[0]=='['){
                status = 2; // for start evt number
            }else{
                status = 3; // for end evt number
            }

            if(!(nameInput[size_]=='{'||nameInput[size_]=='"'||nameInput[size_]=='['
                        ||nameInput[size_]==']'||nameInput[size_]==','||nameInput[size_]==':'||nameInput[size_]=='}'))
                nameOutput += nameInput[size_];
        }
}
