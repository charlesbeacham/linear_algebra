from decimal import Decimal, getcontext

from vector import Vector

getcontext().prec = 30


class Line(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 2

        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()


    def set_basepoint(self):
        try:
            n = self.normal_vector.coordinates
            c = self.constant_term
            basepoint_coords = ['0']*self.dimension

            initial_index = Line.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e


    def __str__(self):

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector.coordinates

        try:
            initial_index = Line.first_nonzero_index(n)
            terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
                     for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output


    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)
        
    def is_parallel_line(self,l2):
        '''
        Determines if two lines are parallel.  They could also be equal.
        Two lines are parallel if their normal vectors are parallel.
        '''
        
        return self.normal_vector.is_parallel(l2.normal_vector)
        
    def __eq__(self,l2):
        '''
        Determines if two lines are equal.  Two lines are equal if the 
        vector connecting a point from each line is parallel to the line.
        i.e. the vector will be orthogonal to both lines' normal vectors.
        
        The lines also will have infinite intersections.
        '''
        
        #this if statement handles the case of the normal vector being the zero vector
        #if both normal vectors are zero vectors, then if the constants are equal the lines
        #are equal.  If not, then they are not equal.
        if self.normal_vector.is_zero():
            if not l2.normal_vector.is_zero():
                return False
            
            else:
                diff = self.constant_term - l2.constant_term
                return MyDecimal(diff).is_near_zero()
        
        elif l2.normal_vector.is_zero():
            return False
        
        if not self.is_parallel_line(l2):
            return False
            
        x0,y0 = self.basepoint, l2.basepoint
        connect_vector = x0.minus(y0)
        
        if connect_vector.is_ortho(self.normal_vector) and connect_vector.is_ortho(l2.normal_vector):
            return True
        else:
            return False 
            
    def get_intersection(self,l2):
        '''
        Find the intersection point of two lines, or return that there is no intersection or infinite intersections.
        
        If lines are not parallel then AD-BC cannot be zero.  Both A & C
        can also not both be zero b/c then both lines would be horizontal 
        and thus parallel.
        '''
        l_1,l_2 = self, l2
        
        if l_1 == l_2:
            print("Both lines are equal.")   
            return
            
        elif l_1.is_parallel_line(l_2):
            print("Both lines are parallel but not equal.")
            return
        
        else:
            if l_1.normal_vector.coordinates[0] == 0:
                l_1,l_2 = l2,self
            
            A,B = l_1.normal_vector.coordinates
            k1 = l_1.constant_term
            C,D = l_2.normal_vector.coordinates
            k2 = l_2.constant_term
            
            x = (D*k1 - B*k2) / (A*D - B*C)
            y = (-C*k1 + A*k2) / (A*D - B*C)
            
            return Vector([x,y])
            
        


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps