#include <KernelDensityEstimation>

#include <random>

int testInit()
{
    typedef double realScalarType;
    typedef kernelDensityEstimation<realScalarType,GAUSS,DEFAULT> kernelDensityEstimationType;
    typedef kernelDensityEstimationType::realMatrixType realMatrixType;

    //random number generator
    std::mt19937 gen(1234);
    std::normal_distribution<realScalarType> normDist;

    //samples and dimensionality
    const size_t numSamples=1000;
    const size_t numDims=2;

    //generate samples
    realMatrixType samples(numSamples,numDims);
    for(size_t i=0;i<numSamples;++i)
    {
        for(size_t j=0;j<numDims;++j)
        {
            samples(i,j) = normDist(gen);
        }
    }

    kernelDensityEstimationType kde(samples);

    return EXIT_SUCCESS;
}

int main(void)
{
    int ret = 0;
    ret += testInit();
    return 0;
}
