% Test cell2vector.m and make_data_vectors.m

vdata = make_data_vectors(0);
assert(isequal(vdata.vecs, {0}))
assert(isequal(vdata.types, {'double'}))
assert(isequal(vdata.dims, {[1 1]}))

vdata = make_data_vectors(1, 2, 3);
assert(isequal(vdata.vecs, {1, 2, 3}))
assert(isequal(vdata.types, {'double', 'double', 'double'}))
assert(isequal(vdata.dims, {[1 1], [1 1], [1 1]}))

vdata = make_data_vectors(1, [2; 3], [4 6 8; 5 7 9]);
assert(isequal(vdata.vecs, {1, [2 3], [4 5 6 7 8 9]}))
assert(isequal(vdata.types, {'double', 'double', 'double'}))
assert(isequal(vdata.dims, {[1 1], [2 1], [2 3]}))

x1 = {1 3 5; 2 4 6};
vdata_x1 = make_data_vectors(x1);
assert(isequal(vdata_x1.vecs, {1:6}))
assert(isequal(vdata_x1.types, {{'double', 'double', 'double'; 'double', 'double', 'double'}}))
assert(isequal(vdata_x1.dims, {{[1 1], [1 1], [1 1]; [1 1], [1 1], [1 1]}}))

x2 = {[7 8 9], 10, [11; 12], [13 15 17; 14 16 18]};
vdata_x2 = make_data_vectors(x2);
assert(isequal(vdata_x2.vecs, {7:18}))
assert(isequal(vdata_x2.types, {{'double', 'double', 'double', 'double'}}))
assert(isequal(vdata_x2.dims, {{[1 3], [1 1], [2 1], [2 3]}}))

vdata_y = make_data_vectors(x1, x2, [19; 20; 21]);
assert(isequal(cell2mat(vdata_y.vecs), 1:21))
assert(isequal(vdata_y.types, {vdata_x1.types{1}, vdata_x2.types{1}, 'double'}))
assert(isequal(vdata_y.dims, {vdata_x1.dims{1}, vdata_x2.dims{1}, [3 1]}))

A = [1 3 5; 2 4 6];
vdata = make_data_vectors(A, x2, 19, [20; 21]);
assert(isequal(cell2mat(vdata.vecs), 1:21))
assert(isequal(vdata.types, {'double', vdata_x2.types{1}, 'double', 'double'}))
assert(isequal(vdata.dims, {size(A), vdata_x2.dims{1}, [1 1], [2 1]}))
