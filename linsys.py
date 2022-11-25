from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):
        '''
        Takes interger values row1 and row2 and swaps those rows by index where row 0 = index 0.
        '''
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
        '''
        Takes coefficient value and a row id and multiplies that row's plane equation through by the coefficient.
        '''        
        n = self[row].normal_vector
        k = self[row].constant_term
        
        new_normal_vector = n.times_scalar(coefficient)
        new_constant_term = k * coefficient
        
        self[row] = Plane(normal_vector=new_normal_vector,constant_term=new_constant_term)
        
    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
               
        n_add, k_add = self[row_to_add].normal_vector, self[row_to_add].constant_term
        n_add_to, k_add_to = self[row_to_be_added_to].normal_vector, self[row_to_be_added_to].constant_term
        
        new_normal_vector = n_add.times_scalar(coefficient).plus(n_add_to)
        new_constant_term = (k_add * coefficient) + k_add_to
        
        self[row_to_be_added_to] = Plane(normal_vector=new_normal_vector,constant_term=new_constant_term) 

    def compute_triangular_form(self):
        '''
        Initial step of computing the triangular form.
        The algorithm uses a for-while loop to loop over each plane.  If the planes first non_zero index is correct algorithm
        will leave that plane alone.-  Otherwise it will first try to swap it for a correct one.  If the correct plane does not exist it will loop over that plane with a while loop
        until it is either correct or it's normal vector is reduced to the zero vector.
        '''
        system = deepcopy(self)
        sys_ind = system.indices_of_first_nonzero_terms_in_each_row()
        dim = system.dimension
  
        for i in range(len(sys_ind)):
            check_var = i if i < dim else -1            
            while (sys_ind[i] != check_var):  
                if sys_ind[i] == -1: #if normal vector is zero first try swapping it with a later plane to eliminate any zero division errors in the for loop below
                    try:
                        b = [x != -1 for x in sys_ind] 
                        index_to_switch = b.index(True,i+1) #find the first Plane after i where the index != -1
                        system.swap_rows(i,index_to_switch)
                        sys_ind = system.indices_of_first_nonzero_terms_in_each_row() #update while loop exit criteria
                    except:
                        break
                
                try:
                    index_to_switch = sys_ind.index(i,i+1) #check if any other planes have the necessary non-zero index for this position
                    system.swap_rows(i,index_to_switch)  #if so switch them
                    sys_ind = system.indices_of_first_nonzero_terms_in_each_row() #update while loop exit criteria 
                except: #if not need to multiply another row and add until the non_zero_ind == i                       
                    for k in range(i):
                        multiple = system[i].normal_vector.coordinates[k]/system[k].normal_vector.coordinates[k]*-1
                        system.add_multiple_times_row_to_row(multiple,k,i)
                        sys_ind = system.indices_of_first_nonzero_terms_in_each_row() #update while loop exit criteria
                        if sys_ind[i] == check_var:
                            break
                                
        return system
        
    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector.coordinates)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices


    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


p0 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p1 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p2 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p3 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')

s = LinearSystem([p0,p1,p2,p3])

print(s.indices_of_first_nonzero_terms_in_each_row())
print('{},{},{},{}'.format(s[0],s[1],s[2],s[3]))
print(len(s))
print(s)

s[0] = p1
print(s)

print(MyDecimal('1e-9').is_near_zero())
print(MyDecimal('1e-11').is_near_zero())