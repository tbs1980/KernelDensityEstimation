#include <KernelDensityEstimation>

template<typename realScalarType>
int testBandwidth(void)
{
    typedef kde::DiagonalBandwidthMatrix<realScalarType> bandwidthType;
    typedef typename bandwidthType::realMatrixType realMatrixType;
    typedef typename bandwidthType::realVectorType realVectorType;
    typedef typename bandwidthType::indexType indexType;
    typedef typename bandwidthType::bandwidhtMatrixType bandwidhtMatrixType;

    const indexType numDims = 1000;
    const indexType numSamples = 2000;

    realMatrixType data = realMatrixType::Random(numSamples,numDims);

    bandwidthType bw(data);

    realVectorType x = realVectorType::Random(numDims);
    realVectorType xTest = x;

    bandwidhtMatrixType HInvSqrt(numDims);

    for(indexType i=0;i<numDims;++i)
    {
        realScalarType count = (realScalarType) numSamples;
        realScalarType mean = data.col(i).mean();
        realScalarType sqMean = data.col(i).squaredNorm()/count;
        realScalarType sigma = std::sqrt(sqMean - mean*mean);

        HInvSqrt(i) = std::pow( realScalarType(4.)/realScalarType(numDims+2),
            -realScalarType(1)/realScalarType(numDims+4) )
            * std::pow(realScalarType(numSamples),realScalarType(1)/realScalarType(numDims+4) )
            / sigma;
    }

    bw.scale(x);

    xTest = xTest.cwiseProduct(HInvSqrt);

    realScalarType eps = std::numeric_limits<realScalarType>::epsilon();

    for(indexType i=0;i<numDims;++i)
    {
        if(std::abs(x(i) - xTest(i)) > eps)
        {
            std::cout<<"std::abs(x(i) - xTest(i)) = "<<std::abs(x(i) - xTest(i))<<std::endl;
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

int main()
{
    int ret = 0;

    ret = testBandwidth<float>();
    ret = testBandwidth<double>();
    ret = testBandwidth<long double>();

    return ret;
}
