#include <KernelDensityEstimation>

template<typename realScalarType>
void createMultiVariateData(const long numDims)
{
    typedef utils::GaussLogPost<realScalarType> GaussLogPostType;
    typedef typename GaussLogPostType::realVectorType realVectorType;
    typedef typename GaussLogPostType::realMatrixType realMatrixType;
    typedef typename GaussLogPostType::indexType indexType;

    //define the Gaussian posterior
    //const indexType numDims=2;
    realMatrixType Mat=realMatrixType::Random(numDims,numDims);
    realMatrixType sigma=Mat*Mat.transpose();;
    realMatrixType sigmaInv = sigma.inverse();//realMatrixType::Identity(numDims,numDims);

    realVectorType mu = realVectorType::Zero(numDims);


    GaussLogPostType glp(mu,sigmaInv,1234l);

    //number of samples
    const indexType numSamples=3000;

    //generate samples
    std::ofstream outFile;
    if(numDims == 1)
    {
        outFile.open("uni-variate-Gauss.dat",std::ios::trunc);
    }
    else if(numDims ==2 )
    {
        outFile.open("bi-variate-Gauss.dat",std::ios::trunc);
    }
    else
    {
        outFile.open("multi-variate-Gauss.dat",std::ios::trunc);
    }
    outFile<<std::scientific;
    realMatrixType samples(numSamples,numDims);
    for(indexType i=0;i<numSamples;++i)
    {
        realVectorType samp = glp.generate();
        samples.row(i) = samp;

        for(indexType j=0;j<numDims;++j)
        {
            outFile<<std::setprecision(10)<<samp(j)<<",";
        }
        //final column is the log-pdf
        outFile<<glp.compute(samp)<<std::endl;
    }

    outFile.close();
}

template<typename realScalarType>
void GaussAllNeighbours()
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
    realMatrixType Mat=realMatrixType::Random(numDims,numDims);
    realMatrixType sigma=Mat*Mat.transpose();;
    realMatrixType sigmaInv = sigma.inverse();//realMatrixType::Identity(numDims,numDims);

    realVectorType mu = realVectorType::Zero(numDims);


    GaussLogPostType glp(mu,sigmaInv,1234l);

    //number of samples
    const indexType numSamples=10000;

    //generate samples
    realMatrixType samples(numSamples,numDims);
    for(indexType i=0;i<numSamples;++i)
    {
        realVectorType samp = glp.generate();
        samples.row(i) = samp;
    }

    // benchmark
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    // compute kde
    kdeType kde(samples);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "# Elapsed time for building indices: " << elapsed_seconds.count() << "s\n";

    const indexType numTestSamples=1000;
    realMatrixType testSamples(numTestSamples,numDims);
    for(indexType i=0;i<numTestSamples;++i)
    {
        realVectorType samp = glp.generate();
        testSamples.row(i) = samp;
    }

    start = std::chrono::system_clock::now();
    for(indexType i=0;i<numTestSamples;++i)
    {
        realVectorType samp = testSamples.row(i);
        realScalarType logPdf = kde.compute(samp);
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout << "# Elapsed time for "<<numTestSamples<< " calls to compute(): "
        << elapsed_seconds.count() << "s\n";
}

int main()
{
    //createMultiVariateData<double>(2);
    GaussAllNeighbours<double>();

    return 0;
}
