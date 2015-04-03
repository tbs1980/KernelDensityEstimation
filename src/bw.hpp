#ifndef KDE_BANDWIDTH_HPP
#define KDE_BANDWIDTH_HPP

namespace kde
{
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
                realScalarType sqMean = dataMatrix.col(i).squaredNorm()/cnt;
                realScalarType sigma = std::sqrt(sqMean - mean*mean);
                bandwidth(i) = sigma*std::pow(realScalarType(3.)*cnt/realScalarType(4.),
                    realScalarType(-1./5.));
            }

        }
    };

    template<typename _realScalarType>
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
            realScalarType alpha = realScalarType(1.)/( realScalarType(2.)*std::sqrt(M_PI) ) ;
            realScalarType sigma = realScalarType(1.);
            realScalarType n = (realScalarType) data.rows();
            realScalarType q = stiffnessIntegral(w,min,max,data);
            return w - std::pow(n*q*std::pow(sigma,4)/alpha,-realScalarType(1./5.)) ;
        }

        static realScalarType stiffnessIntegral(realScalarType const w,
            realScalarType const min, realScalarType const max,
            realVectorType const & data)
        {
            realScalarType epsilon = 1e-4;
            realScalarType cnt(1);
            realScalarType dx = (max - min)/cnt;
            realScalarType curveMax = curvature(max,w,data);;
            realScalarType curveMin = curvature(max,w,data);
            realScalarType yy = 0.5*(curveMax*curveMax + curveMin*curveMin)*dx;
            realScalarType maxN = (max - min)/std::sqrt(epsilon);

            maxN  = maxN > 2048 ? 2048 : maxN;

            for(indexType n=2; n<= maxN; n*=2)
            {
                dx *= 0.5;
                realScalarType y(0);
                for(indexType i = 1; i <= n-1; i +=2)
                {
                    realScalarType curveVal = curvature(min + i*dx, w, data);
                    curveMin = curveVal*curveVal;
                    y += curveMin;
                }
                yy = 0.5*yy + y*dx;
                if(n > 8 && std::abs(y*dx-0.5*yy) < epsilon*yy)
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
            return ((z*z) - 1.0)*GaussPDF(x,m,s)/(s*s);
        }
    };

    class OptimalBandwidth
    {
    public:
        template<class realMatrixType,class realVectorType>
        static void compute(realMatrixType const& dataMatrix,
            realVectorType & bandwidth)
        {
            typedef typename realMatrixType::Scalar realScalarType;
            typedef typename realMatrixType::Index indexType;

            const indexType maxIter = 25;
            const realScalarType epsilon(1e-03);
            const realScalarType extension(3);

            assert(dataMatrix.cols() == bandwidth.rows());

            classicBandwidth::compute(dataMatrix,bandwidth);

            for(indexType i=0;i<bandwidth.rows();++i)
            {
                realScalarType dataMin = dataMatrix.col(i).minCoeff() - extension*bandwidth(i);
                realScalarType dataMax = dataMatrix.col(i).maxCoeff() + extension*bandwidth(i);
                realScalarType x0 = bandwidth(i);
                realScalarType y0 =
                    OptimalBandwidthEquation<realScalarType>::evaluate(x0,
                        dataMin,dataMax,dataMatrix.col(i));
                realScalarType x = 0.8*x0;
                realScalarType y =
                    OptimalBandwidthEquation<realScalarType>::evaluate(x,
                        dataMin,dataMax,dataMatrix.col(i));
                indexType iter = 0;

                while(iter < maxIter)
                {
                    x -= y*(x0-x)/(y0-y);
                    y = OptimalBandwidthEquation<realScalarType>::evaluate(x,
                            dataMin,dataMax,dataMatrix.col(i));
                    if(std::abs(y) < epsilon*y0)
                    {
                        break;
                    }
                }

                bandwidth(i) = x;
            }
        }
    };


}

#endif //KDE_BANDWIDTH_HPP
