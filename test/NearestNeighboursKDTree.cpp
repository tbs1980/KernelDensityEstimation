#include <KernelDensityEstimation>

#include <random>

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
int testNearestNeighboursKDTRee()
{
    typedef GaussLogPost<realScalarType> GaussLogPostType;
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
