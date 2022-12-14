from decimal import Decimal, getcontext

from vector import Vector

getcontext().prec = 30


class Plane(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 3

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

            initial_index = Plane.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
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
            initial_index = Plane.first_nonzero_index(n)
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
        raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)
        
        
    def is_parallel_plane(self,p2):
        '''
        Determines if two planes are parallel.  They could also be equal.
        Two planes are parallel if their normal vectors are parallel.
        '''
        
        return self.normal_vector.is_parallel(p2.normal_vector)
        
    def __eq__(self,p2):
        '''
        Determines if two planes are equal.  Two planes are equal if the 
        vector connecting a point from each plane is parallel to the plane.
        i.e. the vector will be orthogonal to both planes' normal vectors.
        
        The planes also will have infinite intersections.
        '''
        
        #this if statement handles the case of the normal vector being the zero vector
        #if both normal vectors are zero vectors, then if the constants are equal the planes
        #are equal.  If not, then they are not equal.
        if self.normal_vector.is_zero():
            if not p2.normal_vector.is_zero():
                return False
            
            else:
                diff = self.constant_term - p2.constant_term
                return MyDecimal(diff).is_near_zero()
        
        elif p2.normal_vector.is_zero():
            return False
            
        if not self.is_parallel_plane(p2):
            return False
            
        x0,y0 = self.basepoint, p2.basepoint
        connect_vector = x0.minus(y0)
        
        if connect_vector.is_ortho(self.normal_vector) and connect_vector.is_ortho(p2.normal_vector):
            return True
        else:
            return False
        
    


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps