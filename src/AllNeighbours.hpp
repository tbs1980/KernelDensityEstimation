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
        typedef typename realMatrixType::Index indexType;
        typedef Eigen::Matrix<indexType,Eigen::Dynamic,1> neighbourIndexVectorType;

        explicit AllNeighbours(realMatrixType const & data)
        :mNIVect(data.rows()),mNumDims(data.cols())
        {
            assert(data.rows()>1);
            assert(data.cols()>0);

            for(indexType i=0;i<mNIVect.rows();++i)
            {
                mNIVect(i) = i;
            }
        }

        /**
         * \brief A function that return the indices of (all) neighbours
         * \param x the point for which neighbours are sought
         * \return The indices of (all) neighbours
         */
        inline neighbourIndexVectorType indices(realVectorType const& x) const
        {
            assert(mNumDims == x.rows());
            return mNIVect;
        }

    private:
        neighbourIndexVectorType mNIVect;
        indexType mNumDims;

    };

}//namespace kde


#endif //KDE_ALLNEIGHBOURS_HPP
