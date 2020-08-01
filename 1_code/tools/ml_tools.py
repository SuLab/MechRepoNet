import numpy as np

from scipy.sparse import issparse, csc_matrix, csr_matrix

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MaxAbsScaler
from sklearn.model_selection import cross_validate
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.feature_selection import SelectKBest, chi2, RFE, SelectFromModel

from xgboost import XGBClassifier


def cor_selector(X, y, feature_names, num_feats):
    cor_list = []
    # calculate the correlation with y for each feature
    for i in range(X.shape[1]):
        if issparse(X):
            x = X[:, i].A.reshape(len(y))
        else:
            x = X[:, i]
        cor = np.corrcoef(x, y)[0, 1]
        cor_list.append(cor)
    # replace NaN with 0
    cor_list = [0 if np.isnan(i) else i for i in cor_list]
    # feature name
    cor_feature = np.array(feature_names)[np.argsort(np.abs(cor_list))[-num_feats:].tolist()].tolist()
    # feature selection? 0 for not select, 1 for select
    return [True if i in cor_feature else False for i in feature_names]

def chi2_selector(X, y, num_feats):
    this_selector = SelectKBest(chi2, k=num_feats)
    this_selector.fit(X, y)
    return this_selector.get_support()

def rfe_selector(X, y, num_feats, random_state=None):
    this_selector = RFE(estimator=LogisticRegression(C=.1, solver='liblinear', random_state=random_state),
                        n_features_to_select=num_feats, step=.2, verbose=5)
    this_selector.fit(X, y)
    return this_selector.get_support()

def embeded_lr_selector(X, y, num_feats, random_state=None):
    this_selector = SelectFromModel(LogisticRegression(penalty="l1", solver='liblinear', random_state=random_state),
                                    max_features=num_feats)
    this_selector.fit(X, y)

    return this_selector.get_support()

def embeded_rf_selector(X, y, num_feats, n_jobs, random_state=None):
    rfc = RandomForestClassifier(n_estimators=100, max_depth=50, n_jobs=n_jobs, random_state=random_state)
    this_selector = SelectFromModel(rfc, max_features=num_feats)
    this_selector.fit(X, y)
    return this_selector.get_support()

def embeded_xgb_selector(X, y, num_feats, n_jobs=1, random_state=None):
    # XGBoost takes 0 as default random state
    if random_state is None:
        random_state = 0
    # Paramaters optimized for speed, rather than accuracy (as we have 5 other estimators also providing votes)
    xgbc = XGBClassifier(max_depth=5, n_estimators=200, learning_rate=.16, min_child_weight=1, colsample_bytree=.8,
                         n_jobs=n_jobs, random_state=random_state)
    this_selector = SelectFromModel(xgbc, max_features=num_feats)
    this_selector.fit(X, y)
    return this_selector.get_support()


class FeatureSelector(BaseEstimator, TransformerMixin):

    def __init__(self, num_features=100, min_selections=4, n_jobs=1, feature_names=None, always_keep=None,
                 random_state=None):
        self.num_features = num_features
        self.min_selections = min_selections
        self.n_jobs = n_jobs
        self.feature_names = feature_names
        self.always_keep = always_keep
        self.random_state = random_state

    def fit(self, X, y=None):

        X_norm = MaxAbsScaler().fit_transform(X)
        if issparse(X):
            if type(X) != csc_matrix:
                X = X.tocsc()
            X_norm = X_norm.tocsc()

        print('Running Cor')
        cor_support = cor_selector(X, y, self.feature_names, self.num_features)
        print('Running Chi2')
        chi_support = chi2_selector(X_norm, y, self.num_features)
        print('Running RFE')
        rfe_support = rfe_selector(X_norm, y, self.num_features, self.random_state)
        print('Running LR')
        embeded_lr_support = embeded_lr_selector(X_norm, y, self.num_features, self.random_state)
        print('Running RF')
        embeded_rf_support = embeded_rf_selector(X, y, self.num_features,
                                                 n_jobs=self.n_jobs, random_state=self.random_state)
        print('Running XG')
        embeded_xgb_support = embeded_xgb_selector(X, y, self.num_features,
                                                   n_jobs=self.n_jobs, random_state=self.random_state)

        feature_selection_df = pd.DataFrame({'feature':self.feature_names, 'pearson':cor_support, 'chi_2':chi_support,
                                             'rfe':rfe_support, 'logistics':embeded_lr_support,
                                             'random_forest':embeded_rf_support, 'xgboost':embeded_xgb_support})

        feature_selection_df['total'] = np.sum(feature_selection_df, axis=1)
        self.feature_selection_df_ = feature_selection_df

        keep_features = feature_selection_df.query('total >= {}'.format(self.min_selections))['feature'].tolist()

        # Keep the features that we always want (e.g. domain expertise)
        if self.always_keep is not None:
            keep_features.extend(self.always_keep)

        self.keep_features_ = [f for f in self.feature_names if f in keep_features]

        return self

    def transform(self, X, y=None):

        if issparse(X) and type(X) != csc_matrix:
            X = X.tocsc()
        return X[:, [i for i, f in enumerate(self.feature_names) if f in self.keep_features_]]



class MeanScaledArcsinhTransformer(TransformerMixin):

    def fit(self, X, y=None):
        if issparse(X):
            self.initial_mean_ = X.tocoo().tocsc().mean(axis=0).A[0]
        else:
            self.initial_mean_ = X.mean(axis=0)

        # If input was DataFrame, Converts resultant series to ndarray
        try:
            self.initial_mean_ = self.initial_mean_.values
        except:
            pass

        # If inital mean == 0, likely all values were zero
        # this prevents issues later.
        self.initial_mean_[np.where(self.initial_mean_ == 0.0)] = 1

        return self

    def transform(self, X, y=None):
        if issparse(X):
            return np.arcsinh(X.tocoo().tocsc().multiply(self.initial_mean_**-1)).tocsc()
        return np.arcsinh(X / self.initial_mean_)



def get_model_coefs(model, X, f_names):
    """Helper Function to quickly return the model coefs and correspoding fetaure names"""

    # Ensure we have a numpy array for the features
    if type(X) == pd.DataFrame:
        X = X.values


    # Grab the coeffiencts
    coef = model.coef_
    # Some models return a double dimension array, others only a single
    if len(coef) != len(f_names):
        coef = coef[0]

    # insert the intercept
    coef = np.insert(coef, 0, model.intercept_)
    names = np.insert(f_names, 0, 'intercept')

    # Calculate z-score scaled coefficients based on the features
    if issparse(X):
        if type(X) != csc_matrix:
            X = X.tocoo().tocsc()
        z_intercept = coef[0] + sum(coef[1:] * X.mean(axis=0).A[0])
        z_coef = coef[1:] * sparse_std(X, axis=1)
        z_coef = np.insert(z_coef, 0, z_intercept)
    else:
        z_intercept = coef[0] + sum(coef[1:] * X.mean(axis=0))
        z_coef = coef[1:] * X.std(axis=0)
        z_coef = np.insert(z_coef, 0, z_intercept)

    # Return
    return pd.DataFrame([names, coef, z_coef]).T.rename(columns={0:'feature', 1:'coef', 2:'zcoef'})


def sparse_std(data, axis=1):
    """take the standard deviation of a sparse matrix"""

    def get_vec_std(vec):
        return vec.A.std(ddof=1)

    stds = []

    # ensure the correct matrix type for easy row or column subsetting
    if axis==1 and type(data) != csc_matrix:
        data = data.tocoo().tocsc()
    if axis==0 and type(data) != csr_matrix:
        data = data.tocoo().tocsr()

    # Get the std for each vector along the given axis individually
    for i in range(data.shape[axis]):
        if axis==1:
            stds.append(get_vec_std(data.getcol(i)))
        elif axis==0:
            stds.append(get_vec_std(data.getrow(i)))

    return np.array(stds)
