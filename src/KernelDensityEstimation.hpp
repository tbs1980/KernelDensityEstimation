#ifndef KERNELDENSITYESTIMATION_HPP
#define KERNELDENSITYESTIMATION_HPP

enum GaussBandwidthOptimisation
{
    DEFAULT,
    SECANT,
    BISECTION
};

enum Kernel
{
    GAUSS,
    BOX,
    EPANECHNIKOV
};

enum Density
{
    PDF,
    CDF
};

/**
 * \brief A class for computing Kernel Density Estimation
 *
 * Based on https://github.com/timnugent/kernel-density
 */
template<typename _realScalarType>
class kernelDensityEstimationBase
{
public:
    typedef _realScalarType realScalarType;
    typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,Eigen::Dynamic> realMatrixType;
    typedef Eigen::Matrix<realScalarType,1,Eigen::Dynamic> realVectorType;

    kernelDensityEstimationBase(realMatrixType const& dataMatrix)
    :mDataMatrix(dataMatrix)
    {
        assert(dataMatrix.rows()>0 and dataMatrix.cols()>0);
        mBandwidth = realVectorType::Ones(dataMatrix.cols())*realScalarType(-1);
    }

protected:
    realMatrixType mDataMatrix;
    realVectorType mBandwidth;
};

template<typename _realScalarType,int KernelType,int BandOptType>
class kernelDensityEstimation : public kernelDensityEstimationBase<_realScalarType>
{};


template<typename _realScalarType>
class kernelDensityEstimation<_realScalarType,GAUSS,DEFAULT>
    : public kernelDensityEstimationBase<_realScalarType>
{
public:
    typedef kernelDensityEstimationBase<_realScalarType> KDEBaseType;
    typedef typename KDEBaseType::realScalarType realScalarType;
    typedef typename KDEBaseType::realMatrixType realMatrixType;
    typedef typename KDEBaseType::realVectorType realVectorType;

    kernelDensityEstimation(realMatrixType const& dataMatrix)
    :KDEBaseType(dataMatrix)
    {

    }

    realScalarType pdf(realVectorType const & x) const
    {
        assert(x.rows() == KDEBaseType::mDataMatrix.cols());

        for(size_t i=0;i<KDEBaseType::dataMatrix.rows();++i)
        {
            realScalarType a(1);
            for(size_t j=0;j<x.rows();++j)
            {
                a *= pdfAt(x(j),KDEBaseType::mDataMatrix(i,j),KDEBaseType::mBandwidth(j));
            }
        }
    }

private:

    realScalarType pdfAt(realScalarType const x, realScalarType const mu, realScalarType const sigma) const
    {
        realScalarType z = (x - mu)/sigma;
        return std::exp(-realScalarType(0.5)*z*z)
            /(sigma*std::sqrt( realScalarType(2.)*M_PI));
    }
};

#endif //KERNELDENSITYESTIMATION_HPP
