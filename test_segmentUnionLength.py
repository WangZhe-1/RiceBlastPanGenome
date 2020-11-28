'''
Author: your name
Date: 2020-09-26 15:38:55
LastEditTime: 2020-09-26 15:44:47
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: /Pan_genome_github/test_segmentUnionLength.py
'''
def segmentUnionLength(segment_list):
    points=[]
    for segment in segment_list:
        points.append((segment[0],False))
        points.append((segment[1],True))
    points_sorted = sorted(points, key=lambda x: int(x[0]))
    result=0
    counter=0
    for point_index in range(0,len(points_sorted)):
        if counter:
            result=result+(points_sorted[point_index][0]-points_sorted[point_index-1][0])
        if points_sorted[point_index][1]:
            counter=counter-1
        else:
            counter=counter+1
    return result
segment_list=[(2,5),(4,8),(9,12),(3,15)]
print(segmentUnionLength(segment_list))