from stompy.grid import unstructured_grid
import stompy.grid.quad_laplacian as quads
from stompy.spatial import wkb2shp
import matplotlib.pyplot as plt
import stompy.plot.cmap as scmap
from stompy.spatial import field
from stompy import utils
from matplotlib.tri import LinearTriInterpolator
from scipy.interpolate import griddata

import xarray as xr
import six

turbo=scmap.load_gradient('turbo.cpt')
##
all_points=wkb2shp.shp2geom('../data/cbec/cbec-survey/All_Points.shp',
                            target_srs='EPSG:26910')


xyz=np.array( [ np.array(p) for p in all_points['geom']] )
xyz[:,2]=all_points['Z_m']

## 
six.moves.reload_module(unstructured_grid)
six.moves.reload_module(quads)

gen=unstructured_grid.UnstructuredGrid.read_pickle('cbec-survey-interp-grid8.pkl')

##

# hand-tuned psi_scale per cell 
#psi_scales={0:0.25,
#            1:4.0}

# Try instead setting anisotropy using the i,j inputs
psi_scales={0:1.0,
            1:1.0}

dems=[] # accumulate the results

cell_select=[1]
# cell_select=range(gen.Ncells())

for c in cell_select:
    qgs=quads.QuadsGen(gen,cells=[c])

    # for the lagoon:
    # zoom=(552131, 552509., 4124173., 4124631.)
    # For the north channel
    # zoom=(551871., 553048., 4124313, 4125042)

    # or generically:
    zoom=qgs.g_final.bounds()

    as_built_dem=field.GdalGrid('../data/cbec/to_RCD/asbuilt/merged.tif',
                                geo_bounds=zoom)

    #  Boundary nodes of g get data from the raster.
    #  The rest of AllPoints data takes an IJ from the grid.

    # Pick up samples along the boundary from the DEM
    g=qgs.g_final
    boundaries=g.boundary_linestrings()

    xy_boundary=boundaries[0]

    z_boundary=as_built_dem(xy_boundary)
    xyz_boundary=np.c_[ xy_boundary,z_boundary]

    all_xyz=np.concatenate( [ xyz_boundary,
                              xyz],axis=0 )

    if 0:
        plt.figure(1).clf()
        g.plot_edges(color='k',lw=0.5)

        scat=plt.scatter(all_xyz[:,0],all_xyz[:,1],20,all_xyz[:,2],cmap=turbo)
        img=as_built_dem.plot(cmap=turbo,alpha=0.3)

        img.set_clim(scat.get_clim())

    # Aniso interpolation:

    # Inputs:
    #   grid
    #   phi field defined on grid nodes.
    #   psi field defined on grid nodes.

    # samples in x,y,z
    samp_xy=all_xyz[:,:2]
    samp_z=all_xyz[:,2]

    phi=g.nodes['ij'][:,0]
    psi=g.nodes['ij'][:,1]
    tri=g.mpl_triangulation()

    phi_interp=LinearTriInterpolator(tri,phi)
    psi_interp=LinearTriInterpolator(tri,psi)

    samp_phi= phi_interp(samp_xy[:,0],samp_xy[:,1])
    samp_psi= psi_interp(samp_xy[:,0],samp_xy[:,1])

    tgt_x=np.arange(zoom[0],zoom[1],1.0)
    tgt_y=np.arange(zoom[2],zoom[3],1.0)
    tgt_X,tgt_Y=np.meshgrid(tgt_x,tgt_y)

    tgt_phi=phi_interp(tgt_X,tgt_Y)
    tgt_psi=psi_interp(tgt_X,tgt_Y)

    if 0:
        plt.figure(2).clf()
        fig,axs=plt.subplots(1,2,num=2)
        for ax,scal in zip(axs,[samp_phi,samp_psi]):
            g.plot_edges(color='k',lw=0.5,ax=ax)
            scat=ax.scatter(samp_xy[:,0],samp_xy[:,1],20,scal,cmap=turbo)

    if 0:
        plt.figure(3).clf()
        fig,ax=plt.subplots(num=3)
        ax.scatter(samp_phi,samp_psi,30,samp_z,cmap=turbo)

    # 0.2 seemed streaky in the longitudinal direction.
    # 0.4 allows some weird fills on the south.

    # For the lagoon this was okay:
    # psi_scale=0.25
    # For the N channel, that's fine for the short bit,
    # but opposite what we need for the long bit.
    psi_scale=psi_scales[c]

    valid=np.isfinite(samp_phi*samp_psi*samp_z)

    gridded=griddata( np.c_[ samp_phi, psi_scale*samp_psi][valid],
                      samp_z[valid],
                      np.c_[tgt_phi.ravel(),psi_scale*tgt_psi.ravel()] )
    gridded=gridded.reshape( tgt_X.shape )

    z_fld=field.SimpleGrid( F=gridded,
                            extents=zoom )

    z_fld.smooth_by_convolution(kernel_size=5,iterations=3)

    if 1:
        plt.figure(4).clf()

        img=z_fld.plot(cmap=turbo)
        scat=plt.scatter(all_xyz[:,0],all_xyz[:,1],20,all_xyz[:,2],cmap=turbo)

        scat.set_clim(img.get_clim())
        plt.axis(zoom)

    dems.append(z_fld)

    
#  z_fld.write_gdal('lagoon-interp-1m.tif')

##


# Rethinking the intermediate grid
# Currently the target resolution is also used to define the resolution
# of the psi/phi, and the boundary conditions for psi and phi.
# Changing the matrix for phi is a likely solution to get a proper BC for phi
# that will be consistent with the psi BCs
# The psi/phi grid is primarily for creating a consistent orthgonal coordinate system.


##

gen=unstructured_grid.UnstructuredGrid.read_pickle('cbec-survey-interp-grid8.pkl')

gen.delete_cell(0)
gen.delete_orphan_edges()
gen.delete_orphan_nodes()
gen.renumber()
    
qg=quads.QuadGen(gen,execute=False)

# def create_intermediate_grid(self):

# Need two options:
#   one is to generate a grid that respects ij directly
#   another chooses grid spacing that is roughly isotropic
#   to improve the quality of the FD approximations.
# target grid

IJ=np.zeros_like(gen.nodes['ij']) * np.nan

nom_res=4.0
min_steps=2

# Very similar to fill_ij_interp, but we go straight to
# assigning dIJ
cycles=gen.find_cycles(max_cycle_len=1000)

assert len(cycles)==1,"For now, cannot handle multiple cycles"

# Collect the steps so that we can close the sum at the end
for idx in [0,1]:
    steps=[] # [node a, node b, delta]
    for s in cycles:
        # it's a cycle, so we can roll
        is_fixed=np.nonzero( gen.nodes['ij_fixed'][s,idx] )[0]
        if len(is_fixed):
            s=np.roll(s,-is_fixed[0])
            s=np.r_[s,s[0]] # repeat first node at the end
            # Get the new indices for fixed nodes
            is_fixed=np.nonzero( gen.nodes['ij_fixed'][s,idx] )[0]
        
        dists=utils.dist_along( gen.nodes['x'][s] )

        for a,b in zip( is_fixed[:-1],is_fixed[1:] ):
            d_ab=dists[b]-dists[a]
            dij_ab=gen.nodes['ij'][s[b],idx] - gen.nodes['ij'][s[a],idx]
            if dij_ab==0:
                steps.append( [s[a],s[b],0] )
            else:
                n_steps=max(min_steps,d_ab/nom_res)
                dIJ_ab=int( np.sign(dij_ab) * n_steps )
                steps.append( [s[a],s[b],dIJ_ab] )
        steps=np.array(steps)
        err=steps[:,2].sum()

        stepsizes=np.abs(steps[:,2])
        err_dist=np.round(err*np.cumsum(np.r_[0,stepsizes])/stepsizes.sum())
        err_per_step = np.diff(err_dist)
        steps[:,2] -= err_per_step.astype(np.int32)

    # Now steps properly sum to 0.
    IJ[steps[0,0],idx]=0 # arbitrary starting point
    IJ[steps[:-1,1],idx]=np.cumsum(steps[:-1,2])

gen_IJ=gen.copy()
gen_IJ.nodes['ij']=IJ

if 0:
    # Manually change the pigtail coordinates
    # new_j=-5 # maybe worse?
    new_j=6.0 # 5 better..
    for old,new in [ ( [0,0],[0,new_j] ),
                     ( [4,0],[4,new_j] ) ]:
        n=np.nonzero( (gen_IJ.nodes['ij'][:,0]==old[0]) &
                      (gen_IJ.nodes['ij'][:,1]==old[1]) )[0][0]
        gen_IJ.nodes['ij'][n]=new
           

plt.figure(2).clf()
gen_IJ.plot_edges(color='k',lw=0.5)
def ij_label(i,r):
    s=[]
    if r['ij_fixed'][0]:
        s.append( "i=%.1f"%r['ij'][0] )
    if r['ij_fixed'][1]:
        s.append( "j=%.1f"%r['ij'][1] )
    return " ".join(s)
        
gen_IJ.plot_nodes(labeler=ij_label)
plt.axis('tight')
plt.axis('equal')
plt.axis(zoom)

## 
six.moves.reload_module(quads)

qg2=quads.QuadGen(gen=gen_IJ)

qg2.g_final.delete_node_field('ij')
os.unlink('test-channel.nc')
qg2.g_final.write_ugrid('test-channel.nc',overwrite=True)

##

# Even with the evenly distributed grid, some distortion
# is this a problem with BCs?

zoom=(552515., 552627., 4124387., 4124476)

qg.plot_result()
plt.axis(zoom)

qg.plot_psi_phi()
plt.axis(zoom)

##

qg2.plot_result(num=20)
plt.axis(zoom)

qg2.plot_psi_phi(num=21)
plt.axis(zoom)

##

# So the decoupling between BCs is a problem.  Even worse,
# it's not enough to just couple the BCs. Both I and J need
# to be adjusted to best reflect the boundary geometry.

# Methods of setting the BCs:

#   Include dphi/dx ~ dpsi/dy, dphi/dy ~ -dpsi/dx
#   in the matrix, and include a no-tangential gradient BC
#      for phi.

#   Similar, but only in one direction -- don't need this for
#   the no-flux BCs, only the dirichlet BCs.

#   Include a no-tangential gradient BC for phi, and rows
#   that calculate an integral of dpsi

# All of these hinge on having a no-tangential gradient

# What if both psi and phi had a no tangential gradient
# BC, so nobody has dirichlet.  Need something to set
# the overall scale, right?  Each contiguous boundary
# would have an unconstrained degree of freedom.

# There is just one global scale factor
# the additional DOFs might come from the derivatives,
# evaluated at select locations?

##
six.moves.reload_module(quads)

# Need a small test case

# Gamma-shaped convex polygon
# 6 - 4
# | 2 3
# 0 1
gen=unstructured_grid.UnstructuredGrid(max_sides=150,
                                       extra_node_fields=[('ij',np.float64,2)])

# A starting point that needs some tweaking
if 1:
    n0=gen.add_node(x=[0,0]    ,ij=[0,0])
    n1=gen.add_node(x=[100,0]  ,ij=[3,0])
    n2=gen.add_node(x=[100,205],ij=[3,10])
    n3=gen.add_node(x=[155,190],ij=[6,10])
    n4=gen.add_node(x=[175,265],ij=[6,15])
    n6=gen.add_node(x=[0,300]  ,ij=[0,15])
if 0: # Can I do it manually?
    # Pretty close, but have to adjust both
    # i *and* j. Maybe because I'm fixing the
    # boundary, and an arbitrary boundary
    # is not consistent with an arbitrary
    # set of psi BCs.
    n0=gen.add_node(x=[0,0]    ,ij=[0,0])
    n1=gen.add_node(x=[100,0]  ,ij=[3,0])
    n2=gen.add_node(x=[100,205],ij=[3,10])
    n3=gen.add_node(x=[155,190],ij=[6,10])
    n4=gen.add_node(x=[175,265],ij=[6,13])
    n6=gen.add_node(x=[0,300]  ,ij=[0,13])

gen.add_cell_and_edges(nodes=[n0,n1,n2,n3,n4,n6])

# For this, no need to infer ijs like gen_IJ above.  That's a different step
qg=quads.QuadGen(gen=gen)
qg.plot_result()

qg.plot_psi_phi()

##
from scipy import sparse

# First, adapt the calc_psi_phi() and NodeDiscretization to
# implement zero tangential gradients
self=qg

gtri=self.g_int
self.nd=nd=quads.NodeDiscretization(gtri)

e2c=gtri.edge_to_cells()

# new way that computes both in one go

# check boundaries and determine where Laplacian BCs go
boundary=e2c.min(axis=1)<0
i_dirichlet_nodes={} # for psi
j_dirichlet_nodes={} # for phi

i_tan_groups=[]
j_tan_groups=[]
i_tan_groups_i=[]
j_tan_groups_j=[]

if 0: # Original way -- copy i,j to BCs
    for e in np.nonzero(boundary)[0]:
        n1,n2=gtri.edges['nodes'][e]
        i1=gtri.nodes['ij'][n1,0]
        i2=gtri.nodes['ij'][n2,0]
        if i1==i2:
            i_dirichlet_nodes[n1]=i1
            i_dirichlet_nodes[n2]=i2
        j1=gtri.nodes['ij'][n1,1]
        j2=gtri.nodes['ij'][n2,1]
        if j1==j2:
            # So why does this need to be inverted?
            j_dirichlet_nodes[n1]=-j1
            j_dirichlet_nodes[n2]=-j2
else:
    # Try zero tangential nodes.  Current code will be under-determined
    # possibly even with the derivative constraints
    bcycle=gtri.boundary_cycle()
    n1=bcycle[-1]
    i_grp=None
    j_grp=None
    
    for n2 in bcycle:
        i1=gtri.nodes['ij'][n1,0]
        i2=gtri.nodes['ij'][n2,0]
        j1=gtri.nodes['ij'][n1,1]
        j2=gtri.nodes['ij'][n2,1]
        if i1==i2:
            if i_grp is None:
                i_grp=[n1]
                i_tan_groups.append(i_grp)
                i_tan_groups_i.append(i1)
                j_grp=None
            i_grp.append(n2)
        elif j1==j2:
            if j_grp is None:
                j_grp=[n1]
                j_tan_groups.append(j_grp)
                j_tan_groups_j.append(j1)
                i_grp=None
            j_grp.append(n2)
        else:
            print("Don't know how to deal with non-cartesian edges")
        n1=n2

    # Set the range of psi to [-1,1], and pin some j to 1.0
    low_i=np.argmin(i_tan_groups_i)
    high_i=np.argmax(i_tan_groups_i)
    
    i_dirichlet_nodes[i_tan_groups[low_i][0]]=-1
    i_dirichlet_nodes[i_tan_groups[high_i][0]]=1
    j_dirichlet_nodes[j_tan_groups[1][0]]=1
    
# HERE:
#  6 edges lead to 6 zero rows.
#  Say I impose one scale row --
#   and two centering rows
#  with no more, phi is all zero.
#  

Mblocks=[]
Bblocks=[]
if 1: # PSI
    M_psi_Lap,B_psi_Lap=nd.construct_matrix(op='laplacian',
                                            dirichlet_nodes=i_dirichlet_nodes,
                                            zero_tangential_nodes=i_tan_groups)
    Mblocks.append( [M_psi_Lap,None] )
    Bblocks.append( B_psi_Lap )
if 1: # PHI
    M_phi_Lap,B_phi_Lap=nd.construct_matrix(op='laplacian',
                                            dirichlet_nodes=j_dirichlet_nodes,
                                            zero_tangential_nodes=j_tan_groups)
    Mblocks.append( [None,M_phi_Lap] )
    Bblocks.append( B_phi_Lap )
if 1:
    # PHI-PSI relationship
    # With the no-tangential-gradient BCs, this becomes necessary. It probably
    # isn't necessary in the interior, and including it results in an over-specified
    # problem.
    Mdx,Bdx=nd.construct_matrix(op='dx')
    Mdy,Bdy=nd.construct_matrix(op='dy')
    Mblocks.append( [Mdy,-Mdx] )
    Mblocks.append( [Mdx, Mdy] )
    Bblocks.append( np.zeros(Mdx.shape[1]) )
    Bblocks.append( np.zeros(Mdx.shape[1]) )


bigM=sparse.bmat( Mblocks )
rhs=np.concatenate( Bblocks )

psi_phi,*rest=sparse.linalg.lsqr(bigM,rhs)
self.psi=psi_phi[:gtri.Nnodes()]
self.phi=psi_phi[gtri.Nnodes():]

qg.plot_psi_phi(thinning=0.5)

i_stop,n_iter,res_1norm, res_2norm, a_norm, a_cond, ar_norm, xnorm, calc_var = rest


