from dask_jobqueue import PBSCluster
cluster = PBSCluster(
    job_extra=["-l nodes=1:ppn=10","-l mem=6000MB"],
    header_skip=["select"],
    processes=10
)
cluster.scale(jobs=10)
from dask.distributed import Client
client = Client(cluster)
import time
def slow_increment(x):
    time.sleep(1)
    return (x+1)

from dask.distributed import progress
futures=client.map(slow_increment,range(50000))
progress(futures)
print(cluster.job_script())
