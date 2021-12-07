from collection import defaultdict
once = defaultdict(list)
zeros = defaultdict(list)

original = []
unknown = []

for u in unknown:
  for o in original:
    if o[0] == u[0] and o[1] == u[1]:
      zeros[u].append(o)
    elif o[0] == u[0]:
      # we have mutation at index 1
      once[u].append((2, o))
    elif o[1] == u[1]:
      # we have mutation at index 0
      once[u].append((1, o))

for k in zeros.keys():
  for v in zeros[k]:
    for i in range(2, 11):
      if k[i] != v[i]:
        once[u].append((i + 1, o))

for k in once.keys():
  new_list = []
  for start, v in once[k]:
    correct = True
    for i in range(start, 11):
      if k[i] != v[i]:
        correct = False
        break
    if correct:
      # first is a position of the mutation
      new_list.append((start - 1, v))
  if len(new_list) > 0:
    # update the key-value
    once[k] = new_list
    try:
      assert len(new_list) == 1, 'List is longer than 1'
    except:
      print(k, len(new_list))
  else:
    # remove the key
    once.pop(k)


# once -> mutated_seq: [(position, original_sec)]
# should be 1 element in the list
