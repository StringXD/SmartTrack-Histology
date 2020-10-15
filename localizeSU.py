"""
to localize the neuron in coordinates

"""

import os
import h5py
import re
import csv
import numpy as np
import pandas as pd
import selectivity as zpy

def traverse(path):
    for (basepath, dirs, files) in os.walk(path):
        if "cluster_info.tsv" in files:
            yield basepath


def getMetaInfo():
    site_file = r"H:\NP histology\NP tracks revised.csv"
    infoL = pd.read_csv(site_file).astype(
        {"mouse_id": "str", "implanting_date": "str"}
    )[
        [
            "mouse_id",
            "track_id",
            "implanting_date",
            "side",
            "depth",
        ]
    ]
    return infoL


def get_miceid_date_imecno(path):
    dateA = re.compile("((19|20)\\d\\d\\d\\d)\\D")
    imecNo = re.compile("imec(\\d)")
    miceId = re.compile("[M_\\\\](\\d{2})_")

    dateGrps = dateA.search(path)
    imecGrps = imecNo.search(path)
    miceGrps = miceId.search(path)

    if dateGrps and imecGrps and miceGrps:
        # print(miceGrps.group(1),', ',path)
        return (miceGrps.group(1), dateGrps.group(1), imecGrps.group(1))
    else:
        print("unresolved path meta data: ", path)
        input("press Enter to continue")
        return (None, None, None)


def get_bsid_duration_who(path):
    bs_id = 0
    time_s = 0
    files = os.listdir(path)
    for f in files:
        if f.endswith("ap.meta"):
            with open(os.path.join(path, f), "r") as file:
                for line in file:
                    bs_id_grps = re.match("imDatBsc_sn=(\\d{1,3})", line)
                    if bs_id_grps:
                        bs_id = bs_id_grps.group(1)
                    fs_grps = re.match("fileSizeBytes=(\\d+)", line)
                    if fs_grps:
                        time_s = int(fs_grps.group(1)) / 385 / 2 / 30000

    if bs_id == "350":
        who_did = "ZHA"
    elif bs_id == "142":
        who_did = "HEM"
    else:
        who_did = "UNKNOWN"
        print("Unknown BS id!")
        input("press Enter to continue")

    return (bs_id, time_s, who_did)


def imecNo2side(who_did, date, imecNo, mid):
    if date == "191130" and mid == "49":
        return "L"

    if date == "191101" and mid == "26":
        return "R"

    if who_did == "HEM" and (int(date)) >= 191028:
        if imecNo == "1":
            return "R"
        elif imecNo == "0":
            return "L"
    else:
        if imecNo == "1":
            return "L"
        elif imecNo == "0":
            return "R"

    print("Error parsing imec No")
    return "X"


def getTrackInfo(infoL, mice_id, date, imecNo, who_did):
    tinfo = infoL.loc[
        (infoL["mouse_id"] == mice_id)
        & (infoL["implanting_date"] == "20" + date)
        & (infoL["side"] == imecNo2side(who_did, date, imecNo, mice_id)),
        ["track_id","depth"],
    ]
    if len(tinfo):
        return (tinfo.iat[0,0],tinfo.iat[0,1])
    else:
        return (0,0)
    


if __name__ == "__main__":
    infoL = getMetaInfo()

    paths = []
    mids = []
    tids = []
    depths = []
    for path in traverse(r"F:\recordingDataOld\DataSum"):
        print(path)
        (bs_id, time_s, who_did) = get_bsid_duration_who(path)
        (mice_id, date, imec_no) = get_miceid_date_imecno(path)
        (tid,depth) = getTrackInfo(infoL, mice_id, date, imec_no, who_did)
        print(mice_id)
        print(tid)
        paths.append(path)
        mids.append(mice_id)
        tids.append(tid)
        depths.append(depth)

    os.chdir(r"E:\prJ\neuropixels\histology location analysis")
    f = h5py.File("path2tid.hdf5", "w")
    string_dt = h5py.special_dtype(vlen=str)
    d1 = f.create_dataset("/path",data=np.array(paths,dtype=object),dtype=string_dt)
    d2 = f.create_dataset("/mid", data=np.array(list(map(int,mids))),dtype=int)
    d3 = f.create_dataset("/tid", data=np.array(tids),dtype=int)
    d4 = f.create_dataset("/depth", data=np.array(depths),dtype=int)
    f.close()



