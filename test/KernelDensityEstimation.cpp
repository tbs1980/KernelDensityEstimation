#include <KernelDensityEstimation>

template<typename realScalarType>
int testGaussAllNeighbours()
{
    typedef kde::Gaussian<realScalarType> kernelType;
    typedef kde::DiagonalBandwidthMatrix<realScalarType> bandwidthType;
    typedef kde::AllNeighbours<realScalarType> neighboursType;
    typedef kde::KernelDensityEstimator<kernelType,bandwidthType,neighboursType> kdeType;
    typedef typename kdeType::realVectorType realVectorType;
    typedef typename kdeType::realMatrixType realMatrixType;
    typedef typename kdeType::indexType indexType;
    typedef utils::GaussLogPost<realScalarType> GaussLogPostType;

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

    // compute kde
    kdeType kde(samples);

    // test
    std::ofstream outFile;
    outFile.open("testGaussLogPostNNAll.dat",std::ios::trunc);
    const indexType numTestSamples=1000;

    for(indexType i=0;i<numTestSamples;++i)
    {
        realVectorType samp = glp.generate();

        outFile<<glp.compute(samp)<<","<<kde.compute(samp)<<std::endl;
    }

    outFile.close();


    return EXIT_SUCCESS;
}

template<typename realScalarType>
int testGaussNearestNeighboursKDTree()
{
    typedef kde::Gaussian<realScalarType> kernelType;
    typedef kde::DiagonalBandwidthMatrix<realScalarType> bandwidthType;
    typedef kde::NearestNeighboursKDTree<realScalarType> neighboursType;
    typedef kde::KernelDensityEstimator<kernelType,bandwidthType,neighboursType> kdeType;
    typedef typename kdeType::realVectorType realVectorType;
    typedef typename kdeType::realMatrixType realMatrixType;
    typedef typename kdeType::indexType indexType;
    typedef utils::GaussLogPost<realScalarType> GaussLogPostType;

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

    // compute kde
    kdeType kde(samples);

    // test
    std::ofstream outFile;
    outFile.open("testGaussLogPostNNKDTree.dat",std::ios::trunc);
    const indexType numTestSamples=1000;

    for(indexType i=0;i<numTestSamples;++i)
    {
        realVectorType samp = glp.generate();

        outFile<<glp.compute(samp)<<","<<kde.compute(samp)<<std::endl;
    }

    outFile.close();


    return EXIT_SUCCESS;
}

int main()
{
    int ret = 0;

    ret += testGaussAllNeighbours<float>();
    ret += testGaussAllNeighbours<double>();
    ret += testGaussAllNeighbours<long double>();

    ret += testGaussNearestNeighboursKDTree<float>();
    ret += testGaussNearestNeighboursKDTree<double>();
    ret += testGaussNearestNeighboursKDTree<long double>();

    return ret;
}
