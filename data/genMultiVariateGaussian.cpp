#include <KernelDensityEstimation>

template<typename realScalarType>
void createMultiVariateData(const long numDims)
{
    typedef utils::GaussLogPost<realScalarType> GaussLogPostType;
    typedef typename GaussLogPostType::realVectorType realVectorType;
    typedef typename GaussLogPostType::realMatrixType realMatrixType;
    typedef typename GaussLogPostType::indexType indexType;

    //define the Gaussian posterior
    realMatrixType Mat=realMatrixType::Random(numDims,numDims);
    realMatrixType sigma=Mat*Mat.transpose();;
    realMatrixType sigmaInv = sigma.inverse();//realMatrixType::Identity(numDims,numDims);

    realVectorType mu = realVectorType::Zero(numDims);

    GaussLogPostType glp(mu,sigmaInv,1234l);

    //number of samples
    const indexType numSamples=10000;

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

int main()
{
    createMultiVariateData<double>(2);
    return 0;
}
