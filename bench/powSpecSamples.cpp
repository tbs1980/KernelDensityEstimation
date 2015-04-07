#include <KernelDensityEstimation>


template<typename realScalarType>
void GaussAllNbs(std::vector< std::vector<realScalarType> > const & dataIn)
{
    typedef kde::Gaussian<realScalarType> kernelType;
    typedef kde::DiagonalBandwidthMatrix<realScalarType> bandwidthType;
    typedef kde::AllNeighbours<realScalarType> neighboursType;
    typedef kde::KernelDensityEstimator<kernelType,bandwidthType,neighboursType> kdeType;
    typedef typename kdeType::realVectorType realVectorType;
    typedef typename kdeType::realMatrixType realMatrixType;
    typedef typename kdeType::indexType indexType;

    realMatrixType samples(dataIn.size()/2,dataIn[0].size()-1);

    for(size_t i = 0; i<dataIn.size()/2; ++i)
    {
        for(size_t j=0; j<dataIn[0].size()-1; ++j)
        {
            samples(i,j) = dataIn[i][j];
        }
    }

    // benchmark
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    kdeType kde(samples);

    end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "# Elapsed time for building indices: " << elapsed_seconds.count() << "s\n";

    // test samples
    realMatrixType testSamples(dataIn.size() - dataIn.size()/2,dataIn[0].size()-1);
    realVectorType testLogLiks(dataIn.size() - dataIn.size()/2);
    for(size_t i = 0 ; i< testSamples.rows() ; ++i)
    {
        for(size_t j=0; j<testSamples.cols(); ++j)
        {
            //std::cout<<i<<"\t"<<j<<"\t"<<testSamples.rows()<<"\t"<<testSamples.cols()<<std::endl;
            testSamples(i,j) = dataIn[dataIn.size()/2 + i][j];
        }
        //last column is the log likelihoods
        testLogLiks(i) = dataIn[dataIn.size()/2 + i][dataIn[0].size()-1];
    }

    start = std::chrono::system_clock::now();
    for(indexType i=0;i<testLogLiks.rows();++i)
    {
        realVectorType samp = testSamples.row(i);
        realScalarType logLik =  kde.compute(samp);
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout << "# Elapsed time for "<<testLogLiks.rows()<< " calls to compute(): "
        << elapsed_seconds.count() << "s\n";

    // write output
    std::cout<<"# Writing a scatter plot file benchGaussLogPostNNAll.dat"<<std::endl;
    std::ofstream outFile;
    outFile.open("benchGaussLogPostNNAll.dat",std::ios::trunc);

    for(indexType i=0;i<testLogLiks.rows();++i)
    {
        realVectorType samp = testSamples.row(i);
        outFile<<testLogLiks(i)<<","<<kde.compute(samp)<<std::endl;
    }

    outFile.close();
}

int main(void)
{
    typedef double realScalarType;
    //read data
    std::string fileName("/arxiv/projects/zxcorr/mcmcChainsForLikelihoodTest_2D/2dpost.extract.dat");
    std::vector< std::vector<realScalarType> > dataIn;
    utils::readData<realScalarType>(fileName,',',dataIn);

    GaussAllNbs<realScalarType>(dataIn);


    return 0;
}
