#include <KernelDensityEstimation>

template<typename realScalarType>
int testBandwidth(void)
{
    typedef kde::DiagonalBandwidthMatrix<realScalarType> bandwidthType;
    typedef typename bandwidthType::realMatrixType realMatrixType;
    typedef typename bandwidthType::indexType indexType;

    const indexType numDims = 1000;
    const indexType numSamples = 2000;

    realMatrixType data = realMatrixType::Random(numSamples,numDims);

    bandwidthType bw(data);


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
