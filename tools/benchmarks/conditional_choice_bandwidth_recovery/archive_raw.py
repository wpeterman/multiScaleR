from pathlib import Path
import hashlib,csv,zipfile
root=Path("tools/benchmarks/conditional_choice_bandwidth_recovery/results")
shards=("initial_grid_part1","initial_grid_part2","part1","part2")
names=("assessment.csv","profiles.csv","refined.csv","starts.csv")
paths=[root/shard/name for shard in shards for name in names]
if (root/"boundary-starts").is_dir():
    paths += [root/"boundary-starts"/name for name in names]
with zipfile.ZipFile(root/"raw-profiles-and-predictions.zip","w",
        compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive, (root/"raw-archive-index.csv").open("w",newline="") as index:
    table=csv.writer(index)
    table.writerow(("archive_path","bytes","sha256"))
    for p in paths:
        name=p.relative_to(root).as_posix()
        data=p.read_bytes()
        archive.writestr(name,data)
        table.writerow((name,len(data),hashlib.sha256(data).hexdigest()))
with zipfile.ZipFile(root/"raw-profiles-and-predictions.zip") as archive, (root/"raw-archive-index.csv").open() as index:
    for row in csv.DictReader(index):
        payload=archive.read(row["archive_path"])
        assert len(payload)==int(row["bytes"]) and hashlib.sha256(payload).hexdigest()==row["sha256"]
print(len(paths),"raw files verified")
