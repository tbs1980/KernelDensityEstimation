#include <KernelDensityEstimation>

template<typename realScalarType>
int testIndices()
{
    typedef kde::AllNeighbours<realScalarType> neighboursType;
    typedef typename neighboursType::realMatrixType realMatrixType;
    typedef typename neighboursType::realVectorType realVectorType;
    typedef typename neighboursType::indexType indexType;
    typedef typename neighboursType::neighbourIndexVectorType neighbourIndexVectorType;

    const indexType numDims = 1000;
    const indexType numSamples = 2000;

    realMatrixType data = realMatrixType::Random(numSamples,numDims);

    neighboursType nb(data);

    realVectorType x = realVectorType::Random(numDims);

    neighbourIndexVectorType nbInds = nb.indices(x);

    if(nbInds.size() != numSamples)
    {
        return EXIT_FAILURE;
    }

    for(size_t i=0;i<numSamples;++i)
    {
        if(nbInds[i]  != i)
        {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

int main(void)
{
    int ret = 0;

    ret += testIndices<float>();
    ret += testIndices<double>();
    ret += testIndices<long double>();

    return ret;
}
