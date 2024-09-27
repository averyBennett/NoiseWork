using System;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.Xml;
using System.Windows.Forms;

public class Noise
{
    public static double random1D(double x)
    {
        long seedX = (long)(x + 1);

        long state = seedX ^ (seedX << 21);
        state ^= state << 11;
        state ^= state >> 37;
        state ^= state << 5;

        state ^= state >> 20;
        state ^= (state << 17);
        state ^= state >> (int)(seedX % 63 + 1);

        double randomVal = ((double)state) * 0.00000000001;
        randomVal = randomVal - ((int)randomVal);
        if (randomVal < 0)
        {
            randomVal *= -1;
        }
        return randomVal;
    }
    public static double random2D(double x, double z)
    {
        long seedX = (long)(x + 1);
        long seedZ = (long)(z + 1);

        long state = seedX ^ (seedX << 21) ^ seedZ ^ (seedZ << 19);
        state ^= state << 11;
        state ^= state >> 37;
        state ^= state << 5;

        state ^= state >> 20;
        state ^= (state << 17);
        state ^= state >> (int)(seedZ % 63 + 1);
        state ^= state << (int)(seedX % 63 + 1);

        double randomVal = ((double)state) * 0.00000000001;
        randomVal = randomVal - ((int)randomVal);
        if (randomVal < 0)
        {
            randomVal *= -1;
        }
        return randomVal;
    }
    public static double random3D(double x, double z, double y)
    {
        long seedX = (long)(x + 1);
        long seedZ = (long)(z + 1);
        long seedY = (long)(y + 1);

        long state = seedX ^ (seedX << 21) ^ seedZ ^ (seedZ << 19) ^ seedY ^ (seedY << 23);
        state ^= state << 11;
        state ^= state >> 37;
        state ^= state << 5;

        state ^= state >> 20;
        state ^= (state << 17);
        state ^= state >> (int)(seedZ % 63 + 1);
        state ^= state << (int)(seedX % 63 + 1);
        state ^= state >> (int)(seedY % 63 + 1);

        double randomVal = ((double)state) * 0.00000000001;
        randomVal = randomVal - ((int)randomVal);
        if (randomVal < 0)
        {
            randomVal *= -1;
        }
        return randomVal;
    }
    public static double linesBetween(double x, double z)
    {
        //Mirrors this noise from the negative planes (Looks bad where mirrored so be careful)
        x = Math.Abs(x);
        z = Math.Abs(z);

        //Put coordinates in a cell
        int cellSize = 200;
        double cellX = (int)(x / cellSize);
        double cellZ = (int)(z / cellSize);

        //Compare to current cell focus point
        double xFocus = (random3D(cellX, cellZ, 0) * cellSize) + cellSize * cellX;
        double zFocus = (random3D(cellX, cellZ, 1) * cellSize) + cellSize * cellZ;
        double FocusX2;
        double FocusZ2;
        double retVal;
        double compareVal;

        if (x > xFocus)
        {
            FocusX2 = (random3D(cellX + 1, cellZ, 0) * cellSize) + cellSize * (cellX + 1);
            FocusZ2 = (random3D(cellX + 1, cellZ, 1) * cellSize) + cellSize * cellZ;
            retVal = slopeAlignment(x, z, xFocus, zFocus, FocusX2, FocusZ2);
        } else
        {
            FocusX2 = (random3D(cellX - 1, cellZ, 0) * cellSize) + cellSize * (cellX - 1);
            FocusZ2 = (random3D(cellX - 1, cellZ, 1) * cellSize) + cellSize * cellZ;
            retVal = slopeAlignment(x, z, xFocus, zFocus, FocusX2, FocusZ2);
        }

        if(z > zFocus)
        {
            FocusX2 = (random3D(cellX, cellZ + 1, 0) * cellSize) + cellSize * cellX;
            FocusZ2 = (random3D(cellX, cellZ + 1, 1) * cellSize) + cellSize * (cellZ + 1);
            compareVal = slopeAlignment(x, z, xFocus, zFocus, FocusX2, FocusZ2);
            if (compareVal > retVal)
            {
                retVal = compareVal;
            }
        }
        else
        {
            FocusX2 = (random3D(cellX, cellZ - 1, 0) * cellSize) + cellSize * cellX;
            FocusZ2 = (random3D(cellX, cellZ - 1, 1) * cellSize) + cellSize * (cellZ - 1);
            compareVal = slopeAlignment(x, z, xFocus, zFocus, FocusX2, FocusZ2);
            if (compareVal > retVal)
            {
                retVal = compareVal;
            }
        }

        if (retVal > 1 || retVal < 0)
        {
            Debug.WriteLine("lastCheck-(" + x + ", " + z + ") : " + retVal);
            retVal = 1;
        }

        return retVal;
    }
    //Rewrite slopeAlignment to create less variables, currently too inefficient
    public static double slopeAlignment(double x, double z, double x1, double z1, double x2, double z2)
    {
        //Returns number for distance of x,z point to line between points 1 and 2
        /*if((x1 + z1) < (x2 + z2))
        {
            double xMid = x1;
            double zMid = z1;
            x1 = x2;
            z1 = z2;
            x2 = xMid;
            z2 = zMid;
        }*/

        double xTox1 = x - x1;
        double xTox2 = x2 - x;
        double zToz1 = z - z1;
        double zToz2 = z2 - z;

        if(Math.Abs(xTox1 + zToz1) > Math.Abs(xTox2 + zToz2))
        {
            double midX = xTox1;
            double midZ = zToz1;
            xTox1 = xTox2;
            zToz1 = zToz2;
            xTox2 = midX;
            zToz2 = midZ;
        }

        /*double horiRatio;
        double retRatio;
        if (xTox2 > xTox1)
        {
            horiRatio = xTox1 / xTox2;
            retRatio = horiRatio * zToz2 - zToz1;
        } else
        {
            horiRatio = xTox2 / xTox1;
            retRatio = horiRatio * zToz1 - zToz2;
        }*/

        double horiRatio = xTox1 / xTox2;
        double retRatio = horiRatio * zToz2 - zToz1;

        if (retRatio < 0)
        {
            retRatio *= -1;
        }

        retRatio = 1 / (1 + retRatio);

        return retRatio;
    }
    public static double slopeAlignment2(double x, double z, double x1, double z1, double x2, double z2)
    {
        double xTox1 = x - x1;
        double xTox2 = x2 - x;
        double zToz1 = z - z1;
        double zToz2 = z2 - z;


        return 1;
    }
    public static double voronoi(double x, double z)
    {
        //Mirrors this noise from the negative planes (Looks bad where mirrored so be careful)
        x = Math.Abs(x);
        z = Math.Abs(z);

        //Put coordinates in a cell
        int cellSize = 200;
        double cellX = (int)(x / cellSize);
        double cellZ = (int)(z / cellSize);

        //Run through the 3 by 3 cell area of focus points
        //double focusCellX = cellX;
        //double focusCellZ = cellZ;
        double shortestDist = cellSize * 2;
        double shortestDist2 = shortestDist * 2;

        for (int i = -1; i <= 1; i++)
        {
            for (int j = -1; j <= 1; j++)
            {
                double currCellX = i + cellX;
                double currCellZ = j + cellZ;

                double xDistance = x - (random3D(currCellX, currCellZ, 0) * cellSize + cellSize * currCellX);
                double zDistance = z - (random3D(currCellX, currCellZ, 1) * cellSize + cellSize * currCellZ);
                double distance = Math.Sqrt(xDistance * xDistance + zDistance * zDistance);

                if (distance < shortestDist)
                {
                    shortestDist2 = shortestDist;
                    shortestDist = distance;
                }
                else if (distance < shortestDist2)
                {
                    shortestDist2 = distance;
                }
            }
        }

        //return random2D(focusCellX, focusCellZ); //Returns cells as chunks
        return (shortestDist) / (shortestDist2); //Returns cell boundaries with bounds as bright
    }
    public static double voronoiNodes(double x, double z)
    {
        //Mirrors this noise from the negative planes (Looks bad where mirrored so be careful)
        x = Math.Abs(x);
        z = Math.Abs(z);

        //Put coordinates in a cell
        int cellSize = 300;
        double cellX = (int)(x / cellSize);
        double cellZ = (int)(z / cellSize);

        //Run through the 3 by 3 cell area of focus points
        double shortestDist = cellSize * 2.0;
        double shortestDist2 = cellSize * 2.0;
        double shortestDist3 = cellSize * 2.0;

        for (int i = -1; i <= 1; i++)
        {
            for (int j = -1; j <= 1; j++)
            {
                double currCellX = i + cellX;
                double currCellZ = j + cellZ;

                double xDistance = x - (random3D(currCellX, currCellZ, 0) * cellSize + cellSize * currCellX);
                double zDistance = z - (random3D(currCellX, currCellZ, 1) * cellSize + cellSize * currCellZ);
                double distance = Math.Sqrt(xDistance * xDistance + zDistance * zDistance);

                if (distance < shortestDist)
                {
                    shortestDist3 = shortestDist2;
                    shortestDist2 = shortestDist;
                    shortestDist = distance;
                }
                else if (distance < shortestDist2)
                {
                    shortestDist3 = shortestDist2;
                    shortestDist2 = distance;
                }
                else if (distance < shortestDist3)
                {
                    shortestDist3 = distance;
                }
            }
        }
        if (shortestDist3 < 0.1)
        {
            shortestDist = 0.1;
            shortestDist2 = 0.1;
            shortestDist3 = 0.1;
        }
        else if (shortestDist2 < 0.1)
        {
            shortestDist = 0.1;
            shortestDist2 = 0.1;
        }
        else if (shortestDist < 0.1)
        {
            shortestDist = 0.1;
        }

        double ratiod = 1 / shortestDist2 + 1 / shortestDist3 + 1 / shortestDist;
        ratiod = ratiod * 30;
        ratiod = 1 / ratiod;
        if(ratiod > 1)
        {
            return 0.5;
        }
        //ratiod = ratiod * ratiod;
        return ratiod * ratiod * ratiod; //Returns cell boundaries with bounds as bright
    }
    public static double voronoiSecond(double x, double z)
    {
        //Mirrors this noise from the negative planes (Looks bad where mirrored so be careful)
        x = Math.Abs(x);
        z = Math.Abs(z);

        //Put coordinates in a cell
        int cellSize = 80;
        double shortestDist = 160;
        double cellX = (int)(x / cellSize);
        double cellZ = (int)(z / cellSize);

        //Run through the 3 by 3 cell area of focus points
        double shortestDist2 = shortestDist;

        for (int i = -1; i <= 1; i++)
        {
            for (int j = -1; j <= 1; j++)
            {
                double currCellX = i + cellX;
                double currCellZ = j + cellZ;

                double xDistance = x - (random3D(currCellX, currCellZ, 0) * cellSize + cellSize * currCellX);
                double zDistance = z - (random3D(currCellX, currCellZ, 1) * cellSize + cellSize * currCellZ);
                double distance = Math.Sqrt(xDistance * xDistance + zDistance * zDistance);

                if (distance < shortestDist)
                {
                    shortestDist2 = shortestDist;
                    shortestDist = distance;
                }
                else if (distance < shortestDist2)
                {
                    shortestDist2 = distance;
                }
            }
        }

        double maxDist = Math.Sqrt(2) * cellSize;
        return shortestDist2 / maxDist; //Returns cell boundaries with bounds as bright
    }
    public static double voronoiMountain(double x, double z)
    {
        //Start with normal voronoiSecond
        //Layer with second smaller voronoiSecond, similar to fractal brownian, but weight based on slope of both layers
        //Idea is to rotate generated pyramids/spikes to have base flush with sides of pyramids on the larger base level...
        //to have cascading lines like the ridges of a mountain slowly going down

        //Get height values and slopes
        double baseVoronoi = voronoiSecond(x, z);
        double slopeLarge = (baseVoronoi - voronoiSecond(x - 0.5, z)) * 2; //Slope over difference of 0.5
        double addedVoronoi = voronoiSecond(x * 2, z * 2) * 0.5;
        double slopeSmall = (addedVoronoi - voronoiSecond(2 * x - 0.5, 2 * z)) * 2; //Slope over difference of 0.5

        //weight added voronoi value
        double weight = Math.Abs(slopeLarge - slopeSmall) / 2;
        if(weight > 1)
        {
            weight = 1;
        }
        addedVoronoi = addedVoronoi * weight;

        int pathOption = 1;
        if(pathOption == 0)
        {
            baseVoronoi += addedVoronoi * (2 - baseVoronoi);
            //baseVoronoi /= 1.5;
        }
        else if (pathOption == 1)
        {
            baseVoronoi -= addedVoronoi;
        } else
        {
            return voronoiSecond(x, z);
        }

        return baseVoronoi;
    }
    public static double perlin(double x, double z)
    {
        //Mirrors this noise from the negative planes (Sharp flip where mirrored so be careful)
        if (x < 0)
        {
            x = -x;
        }
        if (z < 0)
        {
            z = -z;
        }

        //Put coordinates in a cell
        int cellSize = 40; //Chosen Sampling Size
        double baseX = (int)(x / cellSize);
        double baseZ = (int)(z / cellSize);
        double dx = (x % cellSize) / cellSize;
        double dz = (z % cellSize) / cellSize;

        /*if(dx < 0.03 && dz < 0.03)
        {
            return 1.0;
        }*/

        //Calculate corner vectors
        double[] corners = new double[8]; //Using 1D array for efficiency
        for (int n = 0; n < 4; n++)
        {
            int i = n / 2;
            int j = n % 2;
            double cornerX = baseX + i;
            double cornerZ = baseZ + j;
            double randomVal = random2D(cornerX, cornerZ);
            randomVal = randomVal * 2 - 1;
            double flip = Math.Sqrt(1 - (randomVal * randomVal)); //Make unit vector
            flip *= ((random2D(cornerZ, cornerX) > 0.5) ? -1 : 1); //"Randomly" flip direction of unit vector
            corners[n] = randomVal; //Set Vector
            corners[n + 4] = flip;
        }

        //Calculate Corner values, calculates with dot product formula
        double a = corners[0] * dx + corners[4] * dz;
        double b = corners[2] * (dx - 1) + corners[6] * dz;
        double c = corners[3] * (dx - 1) + corners[7] * (dz - 1);
        double d = corners[1] * dx + corners[5] * (dz - 1);

        //Interpolate between corner values
        double val = interpolation(a, b, c, d, dx, dz);

        return (val + 0.71) * 0.7; //Set 0<val<1 (Rounded, true original range +/-sqrt(1/2)=0.7071)
    }
    private static double interpolation(double a, double b, double c, double d, double tx, double tz)
    {
        tx = ((6 * tx - 15) * tx + 10) * tx * tx * tx;
        tz = ((6 * tz - 15) * tz + 10) * tz * tz * tz;//Polynomial interpolation, always t = [0,1), only adjust steepness
        //6t^5 - 15t^4 + 10t^3 = y; 30x^4 - 60x^3 + 30x^2 = dy
        return tx * (-a + b) + tz * (-a + d) + tx * tz * (a - b - d + c) + a;
        //return a*(1 - tx - tz + txz) + b*(tx - txz) + d*(tz - txz) + c*(txz);
        // a - tx a - tz a + txz a + b tx - b txz + d tz - d txz + c txz
        //tx *(-a+b) + tz *(-a+d) + tx * tz *(a-b-d+c) + a
    }
    public static double perlinDerivZ(double x, double z)
    {
        //Mirrors this noise from the negative planes (Sharp flip where mirrored so be careful)
        if (x < 0)
        {
            x = -x;
        }
        if (z < 0)
        {
            z = -z;
        }

        //Put coordinates in a cell
        int cellSize = 40; //Chosen Sampling Size
        double baseX = (int)(x / cellSize);
        double baseZ = (int)(z / cellSize);
        double dx = (x % cellSize) / cellSize;
        double dz = (z % cellSize) / cellSize;

        //Calculate corner vectors
        double[] corners = new double[8]; //Using 1D array for efficiency
        for (int n = 0; n < 4; n++)
        {
            int i = n / 2;
            int j = n % 2;
            double cornerX = baseX + i;
            double cornerZ = baseZ + j;
            double randomVal = random2D(cornerX, cornerZ);
            randomVal = randomVal * 2 - 1;
            double flip = Math.Sqrt(1 - (randomVal * randomVal)); //Make unit vector
            flip *= ((random2D(cornerZ, cornerX) > 0.5) ? -1 : 1); //"Randomly" flip direction of unit vector
            corners[n] = randomVal; //Set Vector
            corners[n + 4] = flip;
        }

        //Calculate Corner values, calculates with dot product formula
        double a = corners[0] * dx + corners[4] * dz;
        double b = corners[2] * (dx - 1) + corners[6] * dz;
        double c = corners[3] * (dx - 1) + corners[7] * (dz - 1);
        double d = corners[1] * dx + corners[5] * (dz - 1);

        //Interpolate between corner values
        double val = interpolationDerivZ(a, b, c, d, dx, dz);

        /*if(val < 0)//Sets to the magnitude of the derivative, with 0 as 0 slope
        {
            val = val * -1;
        }
        return val * 0.51;*/
        return (val + 2.65) * 0.188; //Set 0<val<1 (Rounded, true original range +/-sqrt(1/2)=0.7071)(0 is negative slope)
    }
    public static double perlinDerivX(double x, double z)
    {
        //Mirrors this noise from the negative planes (Sharp flip where mirrored so be careful)
        if (x < 0)
        {
            x = -x;
        }
        if (z < 0)
        {
            z = -z;
        }

        //Put coordinates in a cell
        int cellSize = 40; //Chosen Sampling Size
        double baseX = (int)(x / cellSize);
        double baseZ = (int)(z / cellSize);
        double dx = (x % cellSize) / cellSize;
        double dz = (z % cellSize) / cellSize;

        //Calculate corner vectors
        double[] corners = new double[8]; //Using 1D array for efficiency
        for (int n = 0; n < 4; n++)
        {
            int i = n / 2;
            int j = n % 2;
            double cornerX = baseX + i;
            double cornerZ = baseZ + j;
            double randomVal = random2D(cornerX, cornerZ); //Unit vector x val
            randomVal = randomVal * 2 - 1;
            double flip = Math.Sqrt(1 - (randomVal * randomVal)); //Make unit vector y val
            flip *= ((random2D(cornerZ, cornerX) > 0.5) ? -1 : 1); //"Randomly" flip direction of unit vector
            corners[n] = randomVal; //Set Vector
            corners[n + 4] = flip;
        }

        //Calculate Corner values, calculates with dot product formula
        double a = corners[0] * dx + corners[4] * dz;
        double b = corners[2] * (dx - 1) + corners[6] * dz;
        double c = corners[3] * (dx - 1) + corners[7] * (dz - 1);
        double d = corners[1] * dx + corners[5] * (dz - 1);
        double dAlph = corners[0] + corners[2] + corners[3] + corners[1];
        /*double f = corners[2];
        double g = corners[3];
        double a = corners[1];*/
        /*
        //tx/tz = dx/dz
        (((6 * tx - 15) * tx + 10) * tx * tx * tx) * ((-a + b) + (((6 * tz - 15) * tz + 10) * tz * tz * tz) * (a - b - d + c))
        + (((6 * tz - 15) * tz + 10) * tz * tz * tz) * (-a + d)
        + a
        --
        a = corners[0] * dx + corners[4] * dz; === h * x + d * z
        b = corners[2] * (dx - 1) + corners[6] * dz; = corners[2] * dx - corners[2] + corners[6] * dz === b * x - b + f * z
        c = corners[3] * (dx - 1) + corners[7] * (dz - 1); = corners[3] * dx - corners[3] + corners[7] * dz - corners[7] === c * x - c + g * z - g
        d = corners[1] * dx + corners[5] * (dz - 1); = corners[1] * dx + corners[5] * dz - corners[5] === a * x + e * z - e

        ------
        h = corners[0], a = corners[1], b = corners[2], c = corners[3], d = corners[4], e = corners[5], f = corners[6], g = corners[7]

        extended
        (((6 * x - 15) * x + 10) * x * x * x) * ((-a + b) + (((6 * z - 15) * z + 10) * z * z * z) * (a - b - d + c))
        + (((6 * z - 15) * z + 10) * z * z * z) * (-a + d)
        + a
        --
        //([a] - [b] - [d] + [c]) = (h * x + d * z - b * x + b - f * z - a * x - e * z + e + c * x - c + g * z - g)
        // === x * (h + c - a - b) + z * (d + g - f - e) + b + e - c - g

        (((6 * x - 15) * x + 10) * x * x * x) * (((b - h) * x + (f - d) * z - b) + (((6 * z - 15) * z + 10) * z * z * z) * (x * (h + c - a - b) + z * (d + g - f - e) + b + e - c - g))
        + (((6 * z - 15) * z + 10) * z * z * z) * ((a - h) * x + (e - d) * z - e)
        + h * x + d * z

        (((6 * x - 15) * x + 10) * x * x * x) * (((b - h) * x + (f - d) * z - b) + (((6 * z - 15) * z + 10) * z * z * z) * (x * (h + c - a - b) + z * (d + g - f - e) + b + e - c - g)) + (((6 * z - 15) * z + 10) * z * z * z) * ((a - h) * x + (e - d) * z - e) + h * x + d * z
        (((6 * x - 15) * x + 10) * x * x * x) * (((b - h) * x + (f - d) * z - b) + (((6 * z - 15) * z + 10) * z * z * z) * (x * (h + c - a - b) + z * (d + g - f - j) + b + j - c - g)) + (((6 * z - 15) * z + 10) * z * z * z) * ((a - h) * x + (j - d) * z - j) + h * x + d * z
         
        */

        //Interpolate between corner values
        //double val = interpolationDerivX(a, b, c, d, dx, dz) * e;

        double val = interpolationDerivX3(a, b, c, d, dx, dz);

        /*if (val < 0)//Sets to the magnitude of the derivative, with 0 as 0 slope
        {
            val = val * -1;
        }
        return val * 0.44;*/
        double retval = val * dAlph;
        return retval;
        //return Math.Abs(val / 0.26);
        //return (val + 3.8) * 0.13; //Set 0<val<1 (Rounded, true original range +/-sqrt(1/2)=0.7071)(0 is negative slope)
    }
    private static double interpolationDerivX(double a, double b, double c, double d, double tx, double tz)
    {
        double retVal = (tz * (6 * tz - 15) + 10) * (a + c - b - d);
        retVal += (b - a);
        retVal *= ((((30 * tx) - 60 * tx) + 30) * tx * tx);
        //double dtx = ((30 * tx - 60) * tx + 30) * tx * tx;
        //tx = ((6 * tx - 15) * tx + 10) * tx * tx * tx;
        //tz = ((6 * tz - 15) * tz + 10) * tz * tz * tz;
        //double dtz = ((30 * tz - 60) * tz + 30) * tz * tz;

        //Polynomial interpolation, always t = [0,1), only adjust steepness
        //6t^5 - 15t^4 + 10t^3 = y; 30x^4 - 60x^3 + 30x^2 = dy
        //tx *(-a+b) + tz *(-a+d) + tx * tz *(a-b-d+c) + a
        //a = (tz * (a - b + c - d)) - a + b;
        return retVal;
        //return a * dtx;
    }
    private static double interpolationDerivX3(double a, double b, double c, double d, double tx, double tz)
    {
        double retVal = (tz * (6 * tz - 15) + 10) * (a + c - b - d);
        retVal += (b - a);
        retVal *= ((((30 * tx) - 60 * tx) + 30) * tx * tx);
        //double dtx = ((30 * tx - 60) * tx + 30) * tx * tx;
        //tx = ((6 * tx - 15) * tx + 10) * tx * tx * tx;
        //tz = ((6 * tz - 15) * tz + 10) * tz * tz * tz;
        //double dtz = ((30 * tz - 60) * tz + 30) * tz * tz;

        //Polynomial interpolation, always t = [0,1), only adjust steepness
        //6t^5 - 15t^4 + 10t^3 = y; 30x^4 - 60x^3 + 30x^2 = dy
        //tx *(-a+b) + tz *(-a+d) + tx * tz *(a-b-d+c) + a
        //a = (tz * (a - b + c - d)) - a + b;
        return retVal;
        //return a * dtx;
    }
    private static double interpolationDerivZ(double a, double b, double c, double d, double tx, double tz)
    {
        tx = ((6 * tx - 15) * tx + 10) * tx * tx * tx;
        //tz = ((6 * tz - 15) * tz + 10) * tz * tz * tz;
        double dtz = ((30 * tz - 60) * tz + 30) * tz * tz;

        //Polynomial interpolation, always t = [0,1), only adjust steepness
        //6t^5 - 15t^4 + 10t^3 = y; 30x^4 - 60x^3 + 30x^2 = dy
        a = tx * (a - b + c - d) - a + d;
        return a * dtz;
    }


    public static double perlinDerivX2(double x, double z)
    {
        // Handle negative coordinates mirroring (same logic as perlin)
        if (x < 0)
        {
            x = -x;
        }
        if (z < 0)
        {
            z = -z;
        }

        // Put coordinates in a cell (same logic as perlin)
        int cellSize = 40;
        double baseX = (int)(x / cellSize);
        double baseZ = (int)(z / cellSize);
        double dx = (x % cellSize) / cellSize;
        double dz = (z % cellSize) / cellSize;

        // Skip if close to cell center for efficiency (same logic as perlin)
        if (dx < 0.03 && dz < 0.03)
        {
            return 0.0; // Derivative is 0 at center
        }

        // Calculate corner vectors (same logic as perlin)
        double[] corners = new double[8];
        for (int n = 0; n < 4; n++)
        {
            int i = n / 2;
            int j = n % 2;
            double cornerX = baseX + i;
            double cornerZ = baseZ + j;
            double randomVal = random2D(cornerX, cornerZ);
            randomVal = randomVal * 2 - 1;
            double flip = Math.Sqrt(1 - (randomVal * randomVal)); // Make unit vector
            flip *= ((random2D(cornerZ, cornerX) > 0.5) ? -1 : 1); // Flipped direction
            corners[n] = randomVal;
            corners[n + 4] = flip;
        }

        // Calculate corner values (derivative of dot product)
        double da_dx = corners[4] * dz; // Derivative of a w.r.t. x
        double db_dx = -corners[6] * dz + corners[4] * dz; // Derivative of b w.r.t. x
        double dc_dx = -corners[7] * dz + corners[6] * dz; // Derivative of c w.r.t. x
        double dd_dx = corners[5] * dz; // Derivative of d w.r.t. x

        // Interpolate derivatives between corner values
        double valDerivX = interpolationDerivX2(da_dx, db_dx, dc_dx, dd_dx, dx, dz);

        return valDerivX * (0.7 * 2); // Adjust scaling based on perlin output range
    }
    private static double interpolationDerivX2(double da_dx, double db_dx, double dc_dx, double dd_dx, double tx, double tz)
    {
        // Derivatives of interpolation polynomial (using chain rule)
        tx = ((30 * tx - 60) * tx + 30) * tx * tx; // d/dx (tx^5...)
        tz = ((30 * tz - 60) * tz + 30) * tz * tz;

        double txz = tx * tz;
        return (da_dx * (1 - tx - tz + txz)) + (db_dx * (tx - txz)) + (dd_dx * (tz - txz)) + (dc_dx * (txz));
    }
}
