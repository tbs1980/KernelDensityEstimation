#ifndef KDE_KERNELDENSITYESTIMATION_HPP
#define KDE_KERNELDENSITYESTIMATION_HPP

namespace kde
{
    template<typename _realScalarType>
    class GaussianKDE
    {
    public:
        typedef _realScalarType realScalarType;
        typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,Eigen::Dynamic> realMatrixType;
        typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,1> realVectorType;
        typedef typename realVectorType::Index indexType;


        GaussianKDE(realMatrixType const& dataMatrix)
        :mDataMatrix(dataMatrix),mBandwidth(mDataMatrix.cols())
        {
            classicBandwidth::compute(mDataMatrix,mBandwidth);
            //OptimalBandwidth::compute(mDataMatrix,mBandwidth);

            /*
            for(indexType i=0;i<mBandwidth.rows();++i)
            {
                std::cout<<"bandwidth "<<i<<" = "<<mBandwidth(i)<<std::endl;
            }
            */
        }

        realScalarType PDF(realVectorType const& x)
        {
            realScalarType d(0);
            for(indexType i=0;i<mDataMatrix.rows();++i)
            {
                realScalarType a(1);
                for(indexType j=0;j<x.rows();++j)
                {
                    a *= GaussPDF(x(j),mDataMatrix(i,j),mBandwidth(j));
                }
                d += a;
            }
            return d/mDataMatrix.rows();
        }

    private:
        realMatrixType mDataMatrix;
        realVectorType mBandwidth;
    };

}

#endif //KDE_KERNELDENSITYESTIMATION_HPP
