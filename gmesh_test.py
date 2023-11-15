import gmsh
import numpy as np

import matplotlib.pyplot as plt
# Initialize gmsh session
gmsh.initialize()

rect = gmsh.model.occ.add_rectangle(0.5, 0, 0, 1.0, 1.0)
gmsh.model.occ.synchronize()
v = gmsh.model.occ.revolve([(2, rect)], 0, 0, 0,  0,
                           1, 0, np.pi, recombine=False)
# vv = gmsh.model.occ.extrude([(2, rect)], 0, 0, 0,  0,
#                     1, 0, np.pi, recombine=True)
surface_loop = gmsh.model.occ.add_surface_loop(
    [rect] + [v_[1] for v_ in v if v_[0] == 2])

gmsh.model.occ.synchronize()
vol = gmsh.model.occ.add_volume([surface_loop])
gmsh.model.occ.remove([(3, 1)])
gmsh.model.occ.remove([(2, i) for i in range(1, 6)])
gmsh.model.occ.synchronize()
mesh = gmsh.model.mesh.generate(3)
gmsh.write("test_3D_1.stl")
gmsh.write("test_3D_1.msh")
# data = np.load('/Users/jeromeguterl/development/gitrm/test_167196_contour.npy',
#                allow_pickle=True).tolist()


# class FluxContour():
#     def __init__(self, contour):
#         self.R = contour[:, 0]/1000.0
#         self.Z = contour[:, 1]/1000.0

#     def __iter__(self):
#         for (R_, Z_) in zip(self.R, self.Z):
#             yield (R_, Z_)


# sep = FluxContour(data['sep'])
# outer = FluxContour(data['outer'])

# # adjust last z
# sep.Z[-1] = outer.Z[-1]
# Z_Dimes = sep.Z[-1]
# P_os = gmsh.model.occ.addPoint(sep.R[-1], 0, sep.Z[-1])
# P_is = gmsh.model.occ.addPoint(sep.R[0], 0, sep.Z[0])
# P_oo = gmsh.model.occ.addPoint(outer.R[-1], 0, outer.Z[-1])
# P_io = gmsh.model.occ.addPoint(outer.R[0], 0, outer.Z[0])

# outer_target = gmsh.model.occ.addLine(P_os, P_oo)
# inner_target = gmsh.model.occ.addLine(P_is, P_io)

# gmsh.model.occ.rotate([(1, inner_target)], 0, 0, 0, 0, 0, 1, np.pi/3)
# gmsh.model.occ.rotate([(1, outer_target)], 0, 0, 0, 0, 0, 1, np.pi/3)

# P_o = [gmsh.model.occ.addPoint(R, 0, Z) for (R, Z) in outer]
# P_s = [gmsh.model.occ.addPoint(R, 0, Z) for (R, Z) in sep]

# # line_o = [gmsh.model.occ.addLine(P_o[i], P_o[i+1]) for i in range(len(P_o)-1)]
# # line_s = [gmsh.model.occ.addLine(P_s[i], P_s[i+1]) for i in range(len(P_s)-1)]

# curve_outer = gmsh.model.occ.add_spline([l for l in P_o])
# curve_sep = gmsh.model.occ.add_spline([l for l in P_s])
# # [(1, l) for l in line_o]
# # [(1, l) for l in line_s]
# gmsh.model.occ.rotate([(1, curve_outer)], 0, 0, 0, 0, 0, 1, np.pi/3)
# gmsh.model.occ.rotate([(1, curve_sep)], 0, 0, 0, 0, 0, 1, np.pi/3)


# cloop = gmsh.model.occ.add_curve_loop(
#     [curve_sep, inner_target, curve_outer, outer_target])
# pl = gmsh.model.occ.add_plane_surface([cloop])

# #
# gmsh.model.occ.remove([(0, t) for t in P_o])
# gmsh.model.occ.remove([(0, t) for t in P_s])
# gmsh.model.occ.remove([(0, P_os), (0, P_is)])
# gmsh.model.occ.remove([(0, P_oo), (0, P_io)])
# gmsh.model.occ.remove([(1, inner_target), (1, outer_target)])
# gmsh.model.occ.remove([(1, curve_sep), (1, curve_outer)])
# v = gmsh.model.occ.revolve([(2, pl)], 0, 0, 0, 0, 0, 1, 2*np.pi/3)

# surface_loop = gmsh.model.occ.add_surface_loop(
#     [v_[1] for v_ in v if v_[0] == 2])
# vol = gmsh.model.occ.add_volume([surface_loop])
# # surface_o = [gmsh.model.occ.revolve(
# #     [(1, l)], 0, 0, 0, 0, 0, 1, 2*np.pi) for l in line_o]
# # surface_s = [gmsh.model.occ.revolve(
# #     [(1, l)], 0, 0, 0, 0, 0, 1, 2*np.pi) for l in line_s]
# # surface_o = gmsh.model.occ.revolve(
# #     [(1, curve_outer)], 0, 0, 0, 0, 0, 1, 2*np.pi)
# # surface_s = gmsh.model.occ.revolve(
# #     [(1, curve_sep)], 0, 0, 0, 0, 0, 1, 2*np.pi)
# # surface_outer_target = gmsh.model.occ.revolve(
# #     [(1, outer_target)], 0, 0, 0, 0, 0, 1, 2*np.pi)
# # surface_inner_target = gmsh.model.occ.revolve(
# #     [(1, inner_target)], 0, 0, 0, 0, 0, 1, 2*np.pi)


# # fig, ax = plt.subplots()
# # ax.plot(w[:,0],w[:,1])
# # ax.plot(o[:,0],o[:,1])
# # ax.plot(i[:,0],i[:,1])
# # ax.plot(s[:,0],s[:,1])


# # *colors* is sequence of rgba tuples.
# # *linestyle* is a string or dash tuple. Legal string values are
# # solid|dashed|dashdot|dotted.  The dash tuple is (offset, onoffseq) where
# # onoffseq is an even length tuple of on and off ink in points.  If linestyle
# # is omitted, 'solid' is used.
# # See `matplotlib.collections.LineCollection` for more information.
# # colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# # # line_segments = LineCollection(segs, linewidths=(0.5, 1, 1.5, 2),
# # #                                colors=colors, linestyle='solid')
# # ax.add_collection(line_segments)
# # ax.set_title('Line collection with masked arrays')
# # plt.show()

# # %%


# # Define DiMES cap surface
# # xDiMES = 1.485
# # yDiMES = 0.0
# # zDiMES = Z_Dimes
# # rDiMES = 0.025

# # # Define small dot on DiMES cap surface
# # xsDot = -0.01 + xDiMES
# # ysDot = 0.0 + yDiMES
# # zsDot = 0.0 + zDiMES
# # rsDot = 0.0005

# # # Define large dot on DiMES cap surface
# # xlDot = 0.0 + xDiMES
# # ylDot = 0.0 + yDiMES
# # zlDot = 0.0 + zDiMES
# # rlDot = 0.005

# # # Create 2D surfaces
# # TagDiMES0 = gmsh.model.occ.addDisk(xDiMES, yDiMES, zDiMES, rDiMES, rDiMES)
# # # surface_outer_target_final = gmsh.model.occ.cut([(surface_outer_target[1])], [
# # #     (2, TagDiMES0)], removeTool=False, removeObject=True)
# # TagsDot = gmsh.model.occ.addDisk(xsDot, ysDot, zsDot, rsDot, rsDot)
# # TaglDot = gmsh.model.occ.addDisk(xlDot, ylDot, zlDot, rlDot, rlDot)
# # TagDiMES = gmsh.model.occ.cut([(2, TagDiMES0)], [(2, TagsDot),
# #                                                  (2, TaglDot)], removeObject=True, removeTool=False)
# # surface_loop = gmsh.model.occ.add_surface_loop([surface_outer_target_final[0][0][1], surface_s[1][1], surface_inner_target[1][1], surface_o[1][1]
# #                                                 ])
# # vol = gmsh.model.occ.add_volume([surface_loop])
# # gmsh.model.addPhysicalGroup(3, [vol], name="simu")
# # Synchronize necessary before mesh setup and generation
gmsh.model.occ.synchronize()
# # %%
# # # Set number of elements on the boundary of each dots and DiMES cap
# # gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
# # gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 50)

# # # Prevent very small elements in small dots
# # gmsh.option.setNumber("Mesh.MeshSizeMin", rsDot/3)

# # # Generate 2D mesh
# # gmsh.model.mesh.setSize(v, 0.001)
# mesh = gmsh.model.mesh.generate(3)

# # Launch the GUI to see the results:
gmsh.fltk.run()
# gmsh.write("test_DiMES9.step")
# # Write mesh into a meshio format

# # %%
# # Close gmsh session
# gmsh.finalize()
