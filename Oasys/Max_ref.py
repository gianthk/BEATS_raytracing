from srxraylib.plot.gol import plot

ang = (in_object_1.get_contents("xoppy_data")[:,0]).tolist()
ref = (in_object_1.get_contents("xoppy_data")[:,1]).tolist()
print(ref.index(max(ref)))
# print(ang[ref.index(max(ref))])

print("Max. Refl.: %f (@ %f x)"%(max(ref), ang[ref.index(max(ref))]))


