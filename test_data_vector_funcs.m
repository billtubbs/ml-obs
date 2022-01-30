% Test make_data_vectors.m and unpack_data_vectors.m

%% Test make_data_vectors

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


%% Test unpack_data_vectors

vdata = make_data_vectors(1);
a = unpack_data_vectors(vdata);
assert(a == 1)

vdata = make_data_vectors(1, 2);
[a, b] = unpack_data_vectors(vdata);
assert(isequal([a b], [1 2]))

vdata = make_data_vectors(1, [2; 3], [4 5]);
[a, b, c] = unpack_data_vectors(vdata);
assert(isequal([a b' c], 1:5))

% Cell arrays
a = 1;
b = [2 3; 4 5];
c = {6, [7; 8], 9};
vdata = make_data_vectors(a, b, c);
[a2, b2, c2] = unpack_data_vectors(vdata);
assert(isequal(a, a2))
assert(isequal(b, b2))
assert(isequal(c, c2))

% Nested cell arrays
a = 1;
b = [2 3; 4 5];
c = {6, {[7 8 9], [10; 11]}, 12};
vdata = make_data_vectors(a, b, c);
[a2, b2, c2] = unpack_data_vectors(vdata);
assert(isequal(a, a2))
assert(isequal(b, b2))
assert(isequal(c, c2))
