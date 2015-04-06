#ifndef KDE_NEAREST_NEIGHBOURS_KDTREE_HPP
#define KDE_NEAREST_NEIGHBOURS_KDTREE_HPP

namespace kde
{
    template<typename _realScalarType>
    class NearestNeighboursKDTree
    {
    public:
        typedef _realScalarType realScalarType;
        typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,1> realVectorType;
        typedef Eigen::Matrix<realScalarType,Eigen::Dynamic,Eigen::Dynamic> realMatrixType;
        typedef typename realMatrixType::Index indexType;
        typedef std::vector<size_t> neighbourIndexVectorType;
        typedef std::vector<realScalarType> distVectType;
        typedef nanoflann::KNNResultSet<realScalarType> resultSetType;

        typedef nanoflann::KDTreeEigenMatrixAdaptor< realMatrixType >  KDTreeType;

        explicit NearestNeighboursKDTree(realMatrixType const & data)
        :mKDTree(data.cols(),data,10),mNumDims(data.cols())
        {
            mKDTree.index->buildIndex();
        }

        inline neighbourIndexVectorType indices(realVectorType const& x) const
        {
            assert(x.rows() == mNumDims);

            const size_t numResults = 3;
            neighbourIndexVectorType   retIndices(numResults);
            distVectType outDistsDqr(numResults);

            resultSetType resultSet(numResults);
            resultSet.init(&retIndices[0], &outDistsDqr[0] );

            mKDTree.index->findNeighbors(resultSet, &x(0), nanoflann::SearchParams(10));

            return retIndices;
        }

    private:
        KDTreeType mKDTree;
        indexType mNumDims;
    };
}


#endif //KDE_NEAREST_NEIGHBOURS_KDTREE_HPP
