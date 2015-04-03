#ifndef KDE_ALLNEIGHBOURS_HPP
#define KDE_ALLNEIGHBOURS_HPP

namespace kde
{
    template<typename _realScalarType>
    class AllNeighbours
    {
    public:
        typedef _realScalarType realScalarType;
        typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,1> realVectorType;
        typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,Eigen::Dynamic> realMatrixType;
        typedef realMatrixType::Index indexType;
        typedef Eigen::Matrix<indexType,Eigen::Dynamic,1> neighbourIndexVectorType;

        explicit AllNeighbours(realMatrixType const & data)
        :mNIVect(data.rows())
        {
            assert(data.rows()>1);
            assert(data.cols()>0);

            for(indexType i=0;i<mNIVect.rows();++i)
            {
                mNI(i) = i;
            }
        }

        /**
         * \brief A function that return the indices of (all) neighbours
         * \param x the point for which neighbours are sought
         * \return The indices of (all) neighbours
         */
        inline neighbourIndexVectorType neighbours(realVectorType const& x) const
        {
            return mNIVect;
        }

    private:
        neighbourIndexVectorType mNIVect;
    };

}//namespace kde


#endif //KDE_ALLNEIGHBOURS_HPP
