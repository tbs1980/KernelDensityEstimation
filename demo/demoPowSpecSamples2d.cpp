#include <KernelDensityEstimation>

template<typename realScalarType>
void GaussAllNbs(std::vector< std::vector<realScalarType> > const & dataIn)
{
    typedef kde::Gaussian<realScalarType> kernelType;
    typedef kde::DiagonalBandwidthMatrix<realScalarType> bandwidthType;
    typedef kde::AllNeighbours<realScalarType> neighboursType;
    typedef kde::KernelDensityEstimator<kernelType,bandwidthType,neighboursType> kdeType;
    typedef typename kdeType::realVectorType realVectorType;
    typedef typename kdeType::realMatrixType realMatrixType;
    typedef typename kdeType::indexType indexType;

    assert(dataIn[0].size() == 3);//will work for a 2-d pdf only!

    realMatrixType samples(dataIn.size()/4,dataIn[0].size()-1);

    for(size_t i = 0; i<dataIn.size()/4; ++i)
    {
        for(size_t j=0; j<dataIn[0].size()-1; ++j)
        {
            samples(i,j) = dataIn[i][j];
        }
    }

    // build kde
    kdeType kde(samples);

    // find the bounds of the 2-d pdf
    const realScalarType xMin = samples.col(0).minCoeff();
    const realScalarType xMax = samples.col(0).maxCoeff();
    const realScalarType yMin = samples.col(1).minCoeff();
    const realScalarType yMax = samples.col(1).maxCoeff();

    // create a plottable data
    const indexType numSteps = 100;
    const realScalarType dx = (xMax - xMin)/(realScalarType) numSteps;
    const realScalarType dy = (yMax - yMin)/(realScalarType) numSteps;

    std::ofstream outFile;
    outFile.open("plot-powSpecSamples2D.dat",std::ios::trunc);
    outFile<<std::scientific;
    for(indexType i=0;i<numSteps;++i)
    {
        realScalarType xi = xMin + i*dx;
        for(indexType j=0;j<numSteps;++j)
        {
            realScalarType yi = yMin + j*dy;
            realVectorType samp(2);
            samp(0) = xi;
            samp(1) = yi;
            //outFile<<std::setprecision(10)<<xi<<","<<yi<<","<<kde.compute(samp)<<std::endl;
            outFile<<std::setprecision(10)<<xi<<","<<yi<<","<<kde.computePDF(samp)<<std::endl;
        }
        outFile<<std::endl;
    }
    outFile.close();

}

int main()
{
    typedef double realScalarType;
    //std::string fileName("../data/powSpecSamples2D.dat");
    std::string fileName("../data/logPowSpecSamples2D.dat");
    std::vector< std::vector<realScalarType> > dataIn;
    utils::readData<realScalarType>(fileName,',',dataIn);

    GaussAllNbs<realScalarType>(dataIn);
}
