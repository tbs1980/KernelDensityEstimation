#include <KernelDensityEstimation>

template<typename realScalarType>
int testGaussian()
{
    typedef kde::Gaussian<realScalarType> kernelType;
    typedef typename kernelType::realVectorType realVectorType;
    typedef typename realVectorType::Index indexType;

    const indexType numDims = 1000;
    realVectorType x = realVectorType::Random(numDims);

    realScalarType sqSum = 0;
    for(indexType i=0;i<x.rows();++i)
    {
        sqSum += x(i)*x(i);
    }

    realScalarType kernelVal = std::exp(-0.5*sqSum)
        *std::pow(realScalarType(2.*M_PI),realScalarType(-numDims*0.5));

    realScalarType kernelValTest = kernelType::compute(x);

    realScalarType eps = std::numeric_limits<realScalarType>::epsilon();

    if(std::abs(kernelVal - kernelValTest) > eps)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int main()
{
    int ret = 0;

    ret += testGaussian<float>();
    ret += testGaussian<double>();
    ret += testGaussian<long double>();

    return 0;
}
