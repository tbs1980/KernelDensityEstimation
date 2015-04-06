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

        typedef nanoflann::RadiusResultSet<realScalarType,size_t> resultSetType;
        typedef nanoflann::KDTreeEigenMatrixAdaptor< realMatrixType >  KDTreeType;

        NearestNeighboursKDTree(realMatrixType const & data)
        :mKDTree(data.cols(),data,10),mNumDims(data.cols()),mRadius(1.)
        {
            mKDTree.index->buildIndex();
            mRadius = getRadius(data);
        }

        neighbourIndexVectorType indices(realVectorType const& x)
        {
            assert(x.rows() == mNumDims);

            return knnSearch(x);
        }

    private:

        neighbourIndexVectorType knnSearch(realVectorType const& x)
        {
            const size_t num_results = 2000;
            neighbourIndexVectorType   ret_indexes(num_results);
            std::vector<realScalarType> out_dists_sqr(num_results);
            nanoflann::KNNResultSet<realScalarType> resultSet(num_results);
            resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
            mKDTree.index->findNeighbors(resultSet, &x(0), nanoflann::SearchParams(10));
            return ret_indexes;
        }

        neighbourIndexVectorType radiusSearch(realVectorType const& x)
        {
            std::vector<std::pair<size_t,realScalarType> > indices_dists;
            resultSetType resultSet(mRadius,indices_dists);

            mKDTree.index->findNeighbors(resultSet, &x(0), nanoflann::SearchParams(10));

            neighbourIndexVectorType retIndices(indices_dists.size());
            for(size_t i=0;i<indices_dists.size();++i)
            {
                retIndices[i] = indices_dists[i].first;
            }
            return retIndices;

        }

        realScalarType getRadius(realMatrixType const & data)
        {
            realVectorType sigma(data.cols());

            for(indexType i=0;i<data.cols();++i)
            {
                //find std-dvn of each dimension
                realScalarType numDims = (realScalarType) data.cols();
                realScalarType count = (realScalarType) data.rows();
                realScalarType mean = data.col(i).mean();
                realScalarType sqMean = data.col(i).squaredNorm()/count;
                sigma(i) = std::sqrt(sqMean - mean*mean);
            }

            realScalarType sigmaVal = sigma.maxCoeff();
            return std::sqrt(sigma.rows())*sigmaVal;
        }

        KDTreeType mKDTree;
        indexType mNumDims;
        realScalarType mRadius;
    };
}


#endif //KDE_NEAREST_NEIGHBOURS_KDTREE_HPP
