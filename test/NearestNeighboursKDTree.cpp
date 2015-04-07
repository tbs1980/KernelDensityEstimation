#include <KernelDensityEstimation>


template<typename realScalarType>
int testNearestNeighboursKDTRee()
{
    typedef utils::GaussLogPost<realScalarType> GaussLogPostType;
    typedef typename GaussLogPostType::realMatrixType realMatrixType;
    typedef typename GaussLogPostType::realVectorType realVectorType;
    typedef typename GaussLogPostType::indexType indexType;
    typedef kde::NearestNeighboursKDTree<realScalarType> neighboursType;
    typedef typename neighboursType::neighbourIndexVectorType neighbourIndexVectorType;

    //define the Gaussian posterior
    const indexType numDims=10;
    realMatrixType sigmaInv = realMatrixType::Identity(numDims,numDims);

    realVectorType mu = realVectorType::Zero(numDims);


    GaussLogPostType glp(mu,sigmaInv,1234l);

    //number of samples
    const indexType numSamples=3000;

    //generate samples
    realMatrixType samples(numSamples,numDims);
    for(indexType i=0;i<numSamples;++i)
    {
        realVectorType samp = glp.generate();
        samples.row(i) = samp;
    }

    neighboursType nbkdt(samples);

    realVectorType query = glp.generate();

    neighbourIndexVectorType nbInds = nbkdt.indices(query);

    /*
    for(size_t i=0;i<nbInds.size();++i)
    {
        std::cout<<i<<"\t"<<nbInds[i]<<std::endl;
    }
    */

    return EXIT_SUCCESS;
}

int main()
{
    int ret = 0;

    ret += testNearestNeighboursKDTRee<float>();
    ret += testNearestNeighboursKDTRee<double>();
    ret += testNearestNeighboursKDTRee<long double>();

    return ret;
}
