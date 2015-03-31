#include <KernelDensityEstimation>

#include <random>

template<typename _realScalarType>
class GaussPost
{
public:
    typedef _realScalarType realScalarType;
    typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,Eigen::Dynamic> realMatrixType;
    typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,1> realVectorType;
    typedef Eigen::LLT<realMatrixType> LLTType;
    typedef typename realVectorType::Index indexType;

    GaussPost(realVectorType const& mu,realMatrixType const& sigmaInv,unsigned long seed)
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

    realScalarType logPost(realVectorType const& x) const
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

int testInit()
{
    typedef double realScalarType;
    typedef kde::GaussianKDE<realScalarType> kernelDensityEstimationType;
    typedef kernelDensityEstimationType::realMatrixType realMatrixType;
    typedef kernelDensityEstimationType::realVectorType realVectorType;
    typedef GaussPost<realScalarType> GaussPostType;
    typedef realMatrixType::Index indexType;

    //define the Gaussian posterior
    const indexType numDims=2;
    realMatrixType sigmaInv = realMatrixType::Identity(numDims,numDims);
    realVectorType mu = realVectorType::Zero(numDims);
    GaussPostType gp(mu,sigmaInv,1234l);

    //samples and dimensionality
    const indexType numSamples=1000;


    //generate samples
    realMatrixType samples(numSamples,numDims);
    for(indexType i=0;i<numSamples;++i)
    {
        realVectorType samp = gp.generate();
        samples.row(i) = samp;
    }

    kernelDensityEstimationType kde(samples);

    for(indexType i=0;i<numSamples;++i)
    {
        realVectorType samp = gp.generate();
        for(indexType j=0;j<samp.rows();++j)
        {
            std::cout<<samp(j)<<",";
        }
        std::cout<<gp.logPost(samp)<<","<<kde.PDF(samp)<<std::endl;
    }

    return EXIT_SUCCESS;
}

int main(void)
{
    int ret = 0;
    ret += testInit();
    return 0;
}
