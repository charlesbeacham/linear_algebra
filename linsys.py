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
            if sys_ind[i] != check_var:  
                if sys_ind[i] == -1: #if normal vector is zero first try swapping it with a later plane to eliminate any zero division errors in the for loop below                    
                    b = [x != -1 for x in sys_ind] 
                    index_to_switch = LinearSystem.find_switch_index(b,True,i) #find the first Plane after i where the index != -1
                    system.swap_rows(i,index_to_switch)
                    sys_ind = system.indices_of_first_nonzero_terms_in_each_row() #update while loop exit criteria
                    if sys_ind[i] == -1: #if after swap still equals -1 then need to quit loop b/c there is nothing left to switch.
                        break
                
                #try swapping for the correct index
                index_to_switch = LinearSystem.find_switch_index(sys_ind,i,i) #check if any other planes have the necessary non-zero index for this position
                system.swap_rows(i,index_to_switch)  #if so switch them
                sys_ind = system.indices_of_first_nonzero_terms_in_each_row() #update while loop exit criteria 
                
                #if the normal_vector value under consideration is 0 (i.e. sys_ind[i] > check_var, then we need to swap it for any subsequant non-zero normal vector value
                if sys_ind[i] > check_var and check_var != -1:              
                    b = [x <= check_var and x != -1 for x in sys_ind] 
                    index_to_switch = LinearSystem.find_switch_index(b,True,i) #find the first Plane after i where the index != -1 or is < check var 
                    system.swap_rows(i,index_to_switch)
                    sys_ind = system.indices_of_first_nonzero_terms_in_each_row() #update while loop exit criteria
                    if sys_ind[i] > check_var and check_var != -1: #if swap wasn't successful then break since there is nothing left to switch
                        break                
                    
                #finally, once we have confirmed that there is a non-zero value in the normal vector and that we weren't able to swap with the correct index, we can clear the other normal vector values until the current index is the leading index.    
                if sys_ind[i] != check_var: #if not need to multiply another row and add until the non_zero_ind == i 
                    try:
                        for k in range(len(sys_ind)):
                            multiple = system[i].normal_vector.coordinates[k]/system[k].normal_vector.coordinates[k]*-1
                            system.add_multiple_times_row_to_row(multiple,k,i)
                            sys_ind = system.indices_of_first_nonzero_terms_in_each_row() #update while loop exit criteria
                            if sys_ind[i] == check_var:
                                break
                    except:
                        break
                        
        return system
    
    @staticmethod
    def find_switch_index(list_to_search,value_to_find,start_point):
        '''
        if it finds the correct index it will return that index, otherwise returns the original index.
        '''        
        try:
            return_index = list_to_search.index(value_to_find,start_point + 1)
        except:
            return_index = start_point

        return return_index
        
    def compute_rref(self):
        tf = self.compute_triangular_form()
        dim = tf.dimension
        sys_ind = tf.indices_of_first_nonzero_terms_in_each_row()
        
        for i in reversed(range(len(tf))):
            if MyDecimal(sum(tf[i].normal_vector.coordinates)).is_near_zero():
                continue            
                     
            #clear above the variable
            j = sys_ind[i]
            for k in range(i):
                multiple = tf[k].normal_vector.coordinates[j]/tf[i].normal_vector.coordinates[j]*-1
                tf.add_multiple_times_row_to_row(multiple,i,k)
                
            #make leading term == 1
            coefficient = Decimal('1') / tf[i].normal_vector.coordinates[j]
            tf.multiply_coefficient_and_row(coefficient,i)
        
        
        return tf
        
    def compute_ge_solution(self):
        '''
        This function will return the result of the unique solution if there is one.  Otherwise it will print if there is no solution or if there are infinite solutions.  No solution if any 0 = k for k non-zero.  Infinite if there is a variable that is not a leading variable in any equation.  Unique solution if each variable is a lead variable and no 0 = k.
        '''
        rref = self.compute_rref()
        
        #test if there is no solutions    
        nv = [p.normal_vector.is_zero() for p in rref] #True if normal_vector is zero vector
        cons = [abs(p.constant_term) < 1e-10 for p in rref] #true if constant term is zero
        no_solution_test = [x == (True,False) for x in zip(nv,cons)] #True if normal_vector is zero and constant term is not zero
        
        #test if there are infinite solutions
        terms = rref.indices_of_first_nonzero_terms_in_each_row()[0:rref.dimension] #an infinite solution exists if there is a variable that is not a leading variable in any equation.
        infinite_solution_test = [x == -1 for x in terms] #True if both normal_vector and constant_term are both 0
        
        if any(no_solution_test): 
            print("There is no solution")
            
        elif any(infinite_solution_test) or len(terms) < rref.dimension:
            print("There are infinite solutions")
            
        else:
            result = [p.constant_term for p in rref][0:rref.dimension]
            return(Vector(result))
            
    def parm_calc(self):
        rref = self.compute_rref()
        
        direction_vectors = rref.extract_direction_vectors_for_parameterization()
        basepoint = rref.extract_basepoint_for_parameterization()
        
        return Parametrization(basepoint, direction_vectors)
        
    def extract_basepoint_for_parameterization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        
        basepoint_coords = [0] * num_variables
        
        for i,p in enumerate(self.planes):
            pivot_var = pivot_indices[i]
            if pivot_var < 0:
                break
            basepoint_coords[pivot_var] = p.constant_term
            
        return Vector(basepoint_coords)        
        
    def extract_direction_vectors_for_parameterization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        free_variable_indices = set(range(num_variables)) - set(pivot_indices)
        
        direction_vectors = []
        
        for free_var in free_variable_indices:
            vector_coords = [0] * num_variables
            vector_coords[free_var] = 1
            for i,p in enumerate(self.planes):
                pivot_var = pivot_indices[i]
                if pivot_var < 0:
                    break
                vector_coords[pivot_var] = -p.normal_vector.coordinates[free_var]
            direction_vectors.append(Vector(vector_coords))
            
        return direction_vectors
        
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
        
class Parametrization(object):
    BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM_MSG = ("The basepoint and direction vectors should all live in the same dimension.")

    def __init__(self, basepoint, direction_vectors):
    
        self.basepoint = basepoint
        self.direction_vectors = direction_vectors
        self.dimension = self.basepoint.dimension
        
        try:
            for v in direction_vectors:
                assert v.dimension == self.dimension
        
        except AssertionError:
            raise Exception(BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM_MSG)

    def __str__(self):

        output = ''
        for coord in range(self.dimension):
            output += 'x_{} = {} '.format(coord + 1,
                                          round(self.basepoint.coordinates[coord], 3))
            for free_var, vector in enumerate(self.direction_vectors):
                output += '+ {} t_{}'.format(round(vector.coordinates[coord], 3),
                                             free_var + 1)
            output += '\n'
        return output           

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