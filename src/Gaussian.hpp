#ifndef KDE_GAUSSIAN_HPP
#define KDE_GAUSSIAN_HPP

namespace kde
{

    template<realScalarType>
    realScalarType GaussPDF(realScalarType const& x,
        realScalarType const& mu, realScalarType const& sigma)
    {
        realScalarType z = (x - mu)/sigma;
        return std::exp(-realScalarType(0.5)*z*z)
            /(sigma*std::sqrt( realScalarType(2.)*M_PI));
    }

    class classicBandwidth
    {
    public:

        template<class realMatrixType,class realVectorType>
        static void compute(realMatrixType const& dataMatrix,
            realVectorType & bandwidth)
        {
            typedef typename realMatrixType::Scalar realScalarType;
            typedef typename realMatrixType::Index indexType;

            assert(dataMatrix.cols() == bandwidth.rows());

            for(indexType i=0;i<bandwidth.rows();++i)
            {
                realScalarType cnt = (realScalarType) dataMatrix.rows();
                realScalarType mean = dataMatrix.col(i).mean();
                realScalarType sqMean = dataMatrix.col(i).squaredNorm();
                realScalarType sigma = std::sqrt(sqMean - mean*mean);
                bandwidth(i) = sigma*std::pow(realScalarType(3.)*cnt/realScalarType(4.),
                    realScalarType(-1./5.));
            }

        }
    };

    template<typename _realSclarType>
    class OptimalBandwidthEquation
    {
    public:

        typedef _realScalarType realScalarType;
        typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,Eigen::Dynamic> realMatrixType;
        typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,1> realVectorType;
        typedef typename realVectorType::Index indexType;

        static realScalarType evaluate(realScalarType const w, realScalarType const min,
            realScalarType const max, realVectorType const & data)
        {
            realScalarType alpha = realScalarType(1.)/(realScalarType(2.)*std::sqrt(M_PI)) ;
            realScalarType sigma = realScalarType(1.);
            realScalarType n = (realScalarType) data.rows();
            realScalarType q = stiffnessIntegral(w,min,max,data);
            return w - std::pow(n*q*std::pow(sigma,4)/alpha,-realScalarType(1./5.)) ;
        }

        static realScalarType stiffnessIntegral(realScalarType const w,
            realScalarType const min, realScalarType const max,
            realVectorType const & data)
        {
            realScalarType epsilon = ;
            realScalarType cnt(1);
            realScalarType dx = (max - min)/cnt;
            realScalarType curveMax = curvature(max,w,data);;
            realScalarType curveMin = curvature(max,w,data);
            realScalarType yy = 0.5*(curveMax*curveMax + curveMin*curveMin)*dx;
            realScalarType maxN = (max - min)/cnt;

            maxN  = maxN > 2048 ? 2048 : maxN;

            for(indexType n=2; n<= maxN; n*=2)
            {
                dx *= 0.5;
                realScalarType y(0);
                for(indexType i = 1; i <= n-1; i +=2)
                {
                    realScalarType curveVal = curvature(min + i*dx, w, data);
                    curveMin = curveVal*curveVal;
                    y += curv_mn;
                }
                yy = 0.5*yy + y*dx;
                if(n > 8 && std::abs(y*dx-0.5*yy) < eps*yy)
                {
                    break;
                }
            }

            return yy;
        }

        static realScalarType curvature(realScalarType const x,
            realScalarType const w, realVectorType const & data)
        {
            realScalarType y(0);
            for(indexType i=0;i<data.rows();++i)
            {
                y += gaussCurvature(x,data(i),w);
            }
            return y/data.rows();
        }

        static realScalarType gaussCurvature(realScalarType const x,
            realScalarType const m, realScalarType const s)
        {
            realScalarType z = (x - m)/s;
            return ((z*z) - 1.0)*GussPDF(x,m,s)/(s*s);
        }


    }

    class OptimalBandwidth
    {
    public:



        template<class realMatrixType,class realVectorType>
        static void compute(realMatrixType const& dataMatrix,
            realVectorType & bandwidth)
        {
            typedef typename realMatrixType::Scalar realScalarType;
            typedef typename realMatrixType::Index indexType;

            assert(dataMatrix.cols() == bandwidth.rows());

            realVectorType defBw(bandwidth);
            classicBandwidth::compute(dataMatrix,defBw);

            for(indexType i=0;i<bandwidth.rows();++i)
            {
                realScalarType x0 = defBw(i);
            }
        }
    }


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

#endif //KDE_GAUSSIAN_HPP
