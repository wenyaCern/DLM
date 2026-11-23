#include "AdaptivePipiKernel.h"

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <stdexcept>
#include <thread>

namespace {
void CheckEdges(const std::vector<double>& edges, const char* name)
{
    if(edges.size() < 2) throw std::runtime_error(std::string("Invalid kernel axis: ") + name);
    for(size_t i=1; i<edges.size(); ++i)
        if(!std::isfinite(edges[i]) || edges[i] <= edges[i-1])
            throw std::runtime_error(std::string("Non-increasing kernel axis: ") + name);
}

unsigned FoldThreadCount()
{
    const char* value=std::getenv("CATS3D_FOLD_THREADS");
    if(!value || !*value) return 1;
    char* end=nullptr;
    const long count=std::strtol(value,&end,10);
    if(*end || count<1 || count>16)
        throw std::runtime_error("CATS3D_FOLD_THREADS must be an integer from 1 to 16");
    return unsigned(count);
}
}

AdaptivePipiKernel::AdaptivePipiKernel(const std::string& fileName)
{
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
    if(!file || file->IsZombie()) throw std::runtime_error("Cannot open adaptive kernel: " + fileName);
    auto* meta = dynamic_cast<TTree*>(file->Get("KernelMetadata"));
    auto* tree = dynamic_cast<TTree*>(file->Get("KernelByQ"));
    if(!meta || !tree || meta->GetEntries()!=1)
        throw std::runtime_error("Adaptive kernel is missing KernelMetadata or KernelByQ");

    unsigned version=0;
    int pdg1=0, pdg2=0, q1q2=0;
    double mass=0;
    std::vector<double>* qe[3]={nullptr,nullptr,nullptr};
    std::vector<double>* re[3]={nullptr,nullptr,nullptr};
    const char* qNames[3]={"qOutEdgesMeV","qSideEdgesMeV","qLongEdgesMeV"};
    const char* rNames[3]={"rOutEdgesFm","rSideEdgesFm","rLongEdgesFm"};
    meta->SetBranchAddress("formatVersion", &version);
    meta->SetBranchAddress("pdg1", &pdg1);
    meta->SetBranchAddress("pdg2", &pdg2);
    meta->SetBranchAddress("q1q2", &q1q2);
    meta->SetBranchAddress("reducedMassMeV", &mass);
    for(unsigned a=0;a<3;++a){
        if(!meta->GetBranch(qNames[a]) || !meta->GetBranch(rNames[a]))
            throw std::runtime_error("Adaptive kernel axis metadata missing");
        meta->SetBranchAddress(qNames[a], &qe[a]);
        meta->SetBranchAddress(rNames[a], &re[a]);
    }
    meta->GetEntry(0);
    if(version!=1 || pdg1!=211 || pdg2!=211 || q1q2!=1 ||
       std::abs(mass-69.785195) > 0.1)
        throw std::runtime_error("Adaptive kernel physics metadata does not match pi+pi+ Coulomb+QS");
    for(unsigned a=0;a<3;++a){
        if(!qe[a] || !re[a]) throw std::runtime_error("Adaptive kernel axis pointer is null");
        qEdges_[a]=*qe[a]; rEdges_[a]=*re[a];
        CheckEdges(qEdges_[a],qNames[a]); CheckEdges(rEdges_[a],rNames[a]);
        if(qEdges_[a].front()!=0 || std::abs(rEdges_[a].front()+rEdges_[a].back())>1e-6)
            throw std::runtime_error("Adaptive kernel requires positive q and symmetric r axes");
        nQ_[a]=unsigned(qEdges_[a].size()-1);
        nR_[a]=unsigned(rEdges_[a].size()-1);
    }
    nR3_=size_t(nR_[0])*nR_[1]*nR_[2];
    const size_t nQ3=size_t(nQ_[0])*nQ_[1]*nQ_[2];
    if(size_t(tree->GetEntries())!=nQ3)
        throw std::runtime_error("Adaptive kernel is incomplete: q-cell count does not match metadata");
    if(nQ3 > std::numeric_limits<size_t>::max()/nR3_)
        throw std::runtime_error("Adaptive kernel dimensions overflow");
    wavefunction2_.resize(nQ3*nR3_);

    unsigned ix=0,iy=0,iz=0;
    std::vector<float>* values=nullptr;
    tree->SetBranchAddress("iqOut",&ix);
    tree->SetBranchAddress("iqSide",&iy);
    tree->SetBranchAddress("iqLong",&iz);
    tree->SetBranchAddress("wavefunction2",&values);
    std::vector<bool> seen(nQ3,false);
    for(Long64_t entry=0;entry<tree->GetEntries();++entry){
        tree->GetEntry(entry);
        if(ix>=nQ_[0] || iy>=nQ_[1] || iz>=nQ_[2] || !values || values->size()!=nR3_)
            throw std::runtime_error("Adaptive kernel contains an invalid q cell");
        const size_t index=QIndex(ix,iy,iz);
        if(seen[index]) throw std::runtime_error("Adaptive kernel contains a duplicate q cell");
        seen[index]=true;
        for(float value:*values)
            if(!std::isfinite(value) || value<0)
                throw std::runtime_error("Adaptive kernel contains a nonfinite or negative wavefunction");
        std::copy(values->begin(),values->end(),wavefunction2_.begin()+index*nR3_);
    }
}

size_t AdaptivePipiKernel::QIndex(unsigned ix,unsigned iy,unsigned iz) const
{
    return (size_t(ix)*nQ_[1]+iy)*nQ_[2]+iz;
}

unsigned AdaptivePipiKernel::NumQ(unsigned axis) const { return nQ_[axis]; }
double AdaptivePipiKernel::QCenter(unsigned axis,unsigned i) const
{
    return 0.5*(qEdges_[axis][i]+qEdges_[axis][i+1]);
}

std::vector<std::pair<unsigned,double>> AdaptivePipiKernel::AxisOverlaps(
    unsigned axis,double lowMeV,double highMeV) const
{
    std::vector<std::pair<unsigned,double>> result;
    const auto& edges=qEdges_[axis];
    for(unsigned i=0;i<nQ_[axis];++i){
        const double positive=std::max(0.0,std::min(highMeV,edges[i+1])-std::max(lowMeV,edges[i]));
        const double negative=std::max(0.0,std::min(highMeV,-edges[i])-std::max(lowMeV,-edges[i+1]));
        if(positive+negative>0) result.emplace_back(i,positive+negative);
    }
    return result;
}

std::vector<double> AdaptivePipiKernel::Fold(double rout,double rside,double rlong) const
{
    const double radius[3]={rout,rside,rlong};
    std::vector<double> axisWeights[3];
    for(unsigned a=0;a<3;++a){
        if(!(radius[a]>0)) throw std::runtime_error("Adaptive kernel requires positive Gaussian radii");
        axisWeights[a].resize(nR_[a]);
        double sum=0;
        for(unsigned i=0;i<nR_[a];++i){
            const double center=0.5*(rEdges_[a][i]+rEdges_[a][i+1]);
            const double width=rEdges_[a][i+1]-rEdges_[a][i];
            const double w=width*std::exp(-center*center/(4*radius[a]*radius[a]));
            axisWeights[a][i]=w; sum+=w;
        }
        for(double& w:axisWeights[a]) w/=sum;
    }
    std::vector<double> source(nR3_);
    for(unsigned ix=0;ix<nR_[0];++ix)
      for(unsigned iy=0;iy<nR_[1];++iy)
        for(unsigned iz=0;iz<nR_[2];++iz)
          source[(size_t(ix)*nR_[1]+iy)*nR_[2]+iz]=
              axisWeights[0][ix]*axisWeights[1][iy]*axisWeights[2][iz];
    const size_t nQ3=size_t(nQ_[0])*nQ_[1]*nQ_[2];
    std::vector<double> result(nQ3);
    const auto foldRange=[&](size_t first,size_t last){
        for(size_t q=first;q<last;++q){
            const float* wf=wavefunction2_.data()+q*nR3_;
            // Every spatial cell keeps its own source weight.
            double sum0=0, sum1=0, sum2=0, sum3=0;
            size_t r=0;
            for(;r+3<nR3_;r+=4){
                sum0+=double(wf[r])*source[r];
                sum1+=double(wf[r+1])*source[r+1];
                sum2+=double(wf[r+2])*source[r+2];
                sum3+=double(wf[r+3])*source[r+3];
            }
            double sum=(sum0+sum1)+(sum2+sum3);
            for(;r<nR3_;++r) sum+=double(wf[r])*source[r];
            result[q]=sum;
        }
    };
    const unsigned workers=std::min<size_t>(FoldThreadCount(),nQ3);
    if(workers==1){
        foldRange(0,nQ3);
    } else {
        const size_t chunk=(nQ3+workers-1)/workers;
        std::vector<std::thread> threads;
        threads.reserve(workers-1);
        for(unsigned worker=1;worker<workers;++worker){
            const size_t first=std::min(nQ3,size_t(worker)*chunk);
            const size_t last=std::min(nQ3,first+chunk);
            threads.emplace_back(foldRange,first,last);
        }
        foldRange(0,std::min(nQ3,chunk));
        for(auto& thread:threads) thread.join();
    }
    return result;
}
