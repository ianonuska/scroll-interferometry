import zarr, numpy as np, json, sys, time
B="https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/PHerc1667/volumes/20251217075048-2.399um-0.2m-78keV-masked.zarr"
z0,z1,y0,y1,x0,x1 = 18408,18664,2000,3024,6352,7376
g=zarr.open_group(B, mode="r"); a=g["0"]; t=time.time()
blk=np.asarray(a[z0:z1, y0:y1, x0:x1]); print("block", blk.shape, blk.dtype, "nonzero frac %.3f"%(blk>0).mean(), "%.0fs"%(time.time()-t), flush=True)
np.save("/workspace/ink3d/block_2p4um.npy", blk)
out=zarr.open_group("/workspace/ink3d/block_2p4um.zarr", mode="w")
arr=out.create_array("0", shape=blk.shape, chunks=(128,128,128), dtype=blk.dtype); arr[:]=blk
out.attrs["multiscales"]=[{"version":"0.4","axes":[{"name":"z","type":"space"},{"name":"y","type":"space"},{"name":"x","type":"space"}],"datasets":[{"path":"0","coordinateTransformations":[{"type":"scale","scale":[1.0,1.0,1.0]}]}]}]
out.attrs["source_bbox"]=[z0,z1,y0,y1,x0,x1]; out.attrs["voxel_um"]=2.399
# 4x downsampled block for run B (mean pooling), stored as 9.6 um
d=blk.reshape(blk.shape[0]//4,4,blk.shape[1]//4,4,blk.shape[2]//4,4).mean((1,3,5)).astype(np.uint8)
outd=zarr.open_group("/workspace/ink3d/block_9p6um.zarr", mode="w"); ad=outd.create_array("0", shape=d.shape, chunks=(64,128,128), dtype=d.dtype); ad[:]=d
outd.attrs["multiscales"]=out.attrs["multiscales"]; outd.attrs["voxel_um"]=9.596; print("downsampled", d.shape, "BLOCK FETCH DONE", flush=True)
