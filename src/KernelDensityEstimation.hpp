#ifndef KDE_KERNELDENSITYESTIMATION_HPP
#define KDE_KERNELDENSITYESTIMATION_HPP

namespace kde
{
    /**
     * \class KernelDensityEstimator
     *
     * \brief A class for estimating kernel density from a set samples
     *
     * \tparam kernelType type of kernel
     * \tparam bandwidthType type of bandwith
     * \tparam neighboursType type of neighbours
     *
     * This class estimates the kernel density of a set of samples using
     * specified kernel, bandwidth and negighbours. The kernel density can be
     * written as \f$ \hat{f}_{\mathbf{H}}(\mathbf{x})
     * = \frac{1}{n} \sum_{i=1}^n K_{\mathbf{H}}(\mathbf{x} - \mathbf{x}_i) \f$.
     *
     */
    template<class kernelType,class bandwidthType,class neighboursType>
    class KernelDensityEstimator
    {
    public:

        typedef typename kernelType::realScalarType realScalarType;
        typedef typename kernelType::realVectorType realVectorType;
        typedef typename bandwidthType::realMatrixType realMatrixType;
        typedef typename bandwidthType::indexType indexType;
        typedef typename neighboursType::neighbourIndexVectorType neighbourIndexVectorType;

        /**
         * \brief A constructor that sets up the class
         * \param data the input data
         */
        explicit KernelDensityEstimator(realMatrixType const & data)
        :mData(data),mBandwidth(data),mNeighbours(data)
        {
        }

        /**
         * \brief A function that returns the PDF at the spcified point
         * \param x the point at which the PDF is sought
         * \return the value of the PDF
         */
        realScalarType compute(realVectorType const& x)
        {
            realScalarType pdf = 0;
            neighbourIndexVectorType nIVect = mNeighbours.neighbours(x);
            for(indexType i=0;i<nIVect.rows();++i)
            {
                indexType ni = nIVect(i)
                pdf += kernelType::compute((x - mData.col(ni)));
            }
            pdf /= (realScalarType)nIVect.rows();
            return pdf;
        }

    private:
        realMatrixType mData; /**< data */
        bandwidthType mBandwidth; /**< bandwidth matrix */
        neighboursType mNeighbours; /**< neighbours */
    }
}//namespace kde

#endif //KDE_KERNELDENSITYESTIMATION_HPP
