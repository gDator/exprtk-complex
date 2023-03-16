#ifndef FRAKTAL_HPP
#define FRAKTAL_HPP

#define EXPRTK_COMPLEX

#include <math.h>
#include <memory>
#include "SFML/Graphics.hpp"
#include <iostream>
#include <vector>
#include <thread>
#include <iostream>
#include "complex.hpp"
#include "exprtk_complex.hpp"
#include "exprtk.hpp"
#include "interpreter_complex.hpp"


typedef ComplexWrapper<double> Complex;
Interpreter<Complex> interpreter;

void pixelValueRecursive(Complex zn, Complex c, int* count)
{
	if(*count <= m_iIteration_limit)
	{
		(*count)++;  
        // std::cout << zn.real << ", " << zn.imag << std::endl;
        interpreter.setValue("z", zn);
        interpreter.setValue("c", c);
		Complex zn1 = interpreter.getValue();   
		if(abs(zn1) <= .m_dConvergence_limit)
		    pixelValueRecursive(zn1, c, count);
		else
			return;
	}
	else
		return;
}

void calculatePixels(sf::RenderWindow &window)
{
    if(!g_calculated)
	{
		//calculation
		Complex c;
		int iterations = 0;
		Complex z0, zn1;
		z0 = g_settings.start_value;
		for (int x = 0; x < window.getSize().x; x++)
		{			//std::cout << "\n";
			for (int y = 0; y < window.getSize().y; y++)
			{
				iterations = 0;
                c.real(normalize(x, 0, window.getSize().x, g_settings.top_left.real(), g_settings.bottom_right.real()));
				c.imag(normalize(y, 0, window.getSize().y, g_settings.top_left.imag(), g_settings.bottom_right.imag()));
                if(pixel_start_value)
                {
                    z0 = c;
                }
				

				pixelValueRecursive(z0, c, &iterations);
				sf::Color col = getColor(iterations);
				drawPixel(x, y, col);
				
			}
		}
		g_calculated = true;
	}
}

#endif /*FRAKTAL_HPP*/