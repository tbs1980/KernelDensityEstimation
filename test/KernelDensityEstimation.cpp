#include <KernelDensityEstimation>

#include <random>
#include <fstream>
#include <iomanip>

template<typename _realScalarType>
class GaussLogPost
{
public:
    typedef _realScalarType realScalarType;
    typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,Eigen::Dynamic> realMatrixType;
    typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,1> realVectorType;
    typedef Eigen::LLT<realMatrixType> LLTType;
    typedef typename realVectorType::Index indexType;

    GaussLogPost(realVectorType const& mu,realMatrixType const& sigmaInv,unsigned long seed)
    :mMu(mu),mSigmaInv(sigmaInv),mGen(seed)
    {
        assert(mMu.rows() == mSigmaInv.rows() );
        assert(mSigmaInv.rows() == mSigmaInv.cols() );

        //find the matrix sigma and its Cholesky
        LLTType lltOfSigmaInv(mSigmaInv.inverse());

        assert(lltOfSigmaInv.info() == Eigen::Success);

        // find the Cholesky
        mChol = lltOfSigmaInv.matrixL();
    }

    realScalarType compute(realVectorType const& x) const
    {
        return -0.5*(mMu-x).transpose()*mSigmaInv*(mMu-x);
    }

    realVectorType generate()
    {
        // generate a vector from N(0,1)
        realVectorType x(mMu.rows());
        for(indexType i=0;i<x.rows();++i)
        {
            x(i) = mNormDist(mGen);
        }

        // rotates the samples
        x = mMu + mChol*x;

        return x;
    }
private:
    realVectorType mMu;
    realMatrixType mSigmaInv;

    std::mt19937 mGen;
    std::normal_distribution<realScalarType> mNormDist;

    realMatrixType mChol;
};

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
    typedef GaussLogPost<realScalarType> GaussLogPostType;

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
    outFile.open("testLogPost.dat",std::ios::trunc);
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

    return ret;
}
