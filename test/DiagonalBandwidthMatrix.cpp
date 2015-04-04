#include <KernelDensityEstimation>

template<typename realScalarType>
int testBandwidth(void)
{
    typedef kde::DiagonalBandwidthMatrix<realScalarType> bandwidthType;
    typedef typename bandwidthType::realMatrixType realMatrixType;
    typedef typename bandwidthType::indexType indexType;
    typedef typename bandwidthType::bandwidhtMatrixType bandwidhtMatrixType;

    const indexType numDims = 1000;
    const indexType numSamples = 2000;

    realMatrixType data = realMatrixType::Random(numSamples,numDims);

    bandwidthType bw(data);

    bandwidhtMatrixType bwMat(numDims);

    realScalarType norm = 1;
    for(indexType i=0;i<numDims;++i)
    {
        realScalarType count = (realScalarType) numSamples;
        realScalarType mean = data.col(i).mean();
        realScalarType sqMean = data.col(i).squaredNorm()/count;
        realScalarType sigma = (sqMean - mean*mean);

        //std::cout<<count<<"\t"<<mean<<"\t"<<sqMean<<"\t"<<sigma<<std::endl;
        //std::cout<<realScalarType(4.)/realScalarType(numDims+2)<<"\t"
        //    <<-realScalarType(1)/realScalarType(numDims+4)<<"\t"
        //    << std::pow( realScalarType(4.)/realScalarType(numDims+2),
        //        -realScalarType(1)/realScalarType(numDims+4) )<<std::endl;
        std::cout<<realScalarType(numSamples)<<"\t"
            <<realScalarType(1)/realScalarType(numDims+4)<<"\t"
            <<std::pow(realScalarType(numSamples),realScalarType(1)/realScalarType(numDims+4) )
            <<std::endl;
        bwMat(i) = std::pow( realScalarType(4.)/realScalarType(numDims+2),
            -realScalarType(1)/realScalarType(numDims+4) )
            * std::pow(realScalarType(numSamples),realScalarType(1)/realScalarType(numDims+4) )
            / sigma;
        norm *= bwMat(i);
    }

    realScalarType eps = std::numeric_limits<realScalarType>::epsilon();

    std::cout<<norm<<"\t"<<bw.norm()<<std::endl;

    if(std::abs(norm - bw.norm()) > eps)
    {
        std::cout<<norm<<"\t"<<bw.norm()<<std::endl;
        return EXIT_FAILURE;
    }


    return EXIT_SUCCESS;
}

int main()
{
    int ret = 0;

    ret = testBandwidth<float>();
    //ret = testBandwidth<double>();
    //ret = testBandwidth<long double>();

    return ret;
}
