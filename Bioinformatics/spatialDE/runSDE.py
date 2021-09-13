import matplotlib.pyplot as plt
import pandas as pd
import NaiveDE
import SpatialDE

meta = pd.read_csv('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/CS6_sections_246.csv', index_col=0)
counts = pd.read_csv('/Users/christopherpenfold/Desktop/spatialDE/AllExp.csv',index_col=0)

plt.scatter(meta['X'][meta['Loc']==246], meta['Y'][meta['Loc']==246], c='k');
plt.axis('equal');
#plt.show()
plt.savefig('246.png')

norm_expr = NaiveDE.stabilize(counts)
norm_expr2 = norm_expr.T
resid_expr = NaiveDE.regress_out(meta, norm_expr, 'np.log(total_counts)').T

sample_resid_expr = resid_expr.sample(n=15000, axis=1, random_state=1)
X = meta[['X', 'Y']]
results = SpatialDE.run(X, sample_resid_expr)
results.to_csv('/Users/christopherpenfold/Desktop/spatialDE/results.csv') # relative position

plt.scatter(meta['X'], meta['Y'], c=norm_expr2['SOX2']);
plt.title('SOX2')
plt.axis('equal')
plt.savefig('246_SOX2.png')
plt.colorbar(ticks=[]);
plt.show()

plt.scatter(meta['X'], meta['Y'], c=norm_expr2['T']);
plt.title('T')
plt.axis('equal')
plt.savefig('246_T.png')
plt.colorbar(ticks=[]);
plt.show()

plt.scatter(meta['X'], meta['Y'], c=norm_expr2['MIXL1']);
plt.title('MIXL1')
plt.axis('equal')
plt.savefig('246_MIXL1.png')
plt.colorbar(ticks=[]);
plt.show()

X = meta[['X', 'Y', 'Z']]
results = SpatialDE.run(X, sample_resid_expr)
results.to_csv('/Users/christopherpenfold/Desktop/spatialDE/results2.csv') # relative position

plt.scatter(meta['X'], meta['Y'], c=norm_expr2['SOX2']);
plt.title('SOX2')
plt.axis('equal')
plt.savefig('246_SOX2_3d.png')
plt.colorbar(ticks=[]);
plt.show()

plt.scatter(meta['X'], meta['Y'], c=norm_expr2['T']);
plt.title('T')
plt.axis('equal')
plt.savefig('246_T_3d.png')
plt.colorbar(ticks=[]);
plt.show()

plt.scatter(meta['X'], meta['Y'], c=norm_expr2['MIXL1']);
plt.title('MIXL1')
plt.axis('equal')
plt.savefig('246_MIXL1_3d.png')
plt.colorbar(ticks=[]);
plt.show()


