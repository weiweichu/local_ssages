#include <stdexcept>
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Basis.h"
#include "schema.h"

namespace SSAGES
{
    BasisFunction* BasisFunction::Build(const Json::Value& json, const std::string& path, uint nbins)
    {
        auto type = json.get("type","none").asString();
        if(type == "Legendre")
            return Legendre::Build(json, path, nbins);
        else if (type == "Chebyshev")
            return Chebyshev::Build(json, path, nbins);
        else if (type == "Fourier")
            return Fourier::Build(json, path, nbins);
        else
            throw std::invalid_argument("Invalid basis set type \"" + type + "\".");
    }

    Chebyshev* Chebyshev::Build(const Json::Value& json,
                              const std::string& path,
                              uint nbin)
    {
        Json::ObjectRequirement validator;
        Json::Value schema;
        Json::Reader reader;
        
        reader.parse(JsonSchema::ChebyshevBasis, schema);
        validator.Parse(schema, path);

        //Validate Inputs
        validator.Validate(json, path);
        if(validator.HasErrors())
            throw BuildException(validator.GetErrors());

        return new Chebyshev(
                json["polynomial_order"].asInt(),
                json["lower_bound"].asDouble(),
                json["upper_bound"].asDouble(),
                nbin
            );
    }

    Fourier* Fourier::Build(const Json::Value& json,
                            const std::string& path,
                            uint nbin)
    {
        Json::ObjectRequirement validator;
        Json::Value schema;
        Json::Reader reader;
        
        reader.parse(JsonSchema::FourierBasis, schema);
        validator.Parse(schema, path);

        //Validate Inputs
        validator.Validate(json, path);
        if(validator.HasErrors())
            throw BuildException(validator.GetErrors());

        return new Fourier(
                json["polynomial_order"].asInt(),
                json["lower_bound"].asDouble(),
                json["upper_bound"].asDouble(),
                nbin
            );
    }
    Legendre* Legendre::Build(const Json::Value& json,
                              const std::string& path,
                              uint nbin)
    {
        Json::ObjectRequirement validator;
        Json::Value schema;
        Json::Reader reader;
        
        reader.parse(JsonSchema::LegendreBasis, schema);
        validator.Parse(schema, path);

        //Validate Inputs
        validator.Validate(json, path);
        if(validator.HasErrors())
            throw BuildException(validator.GetErrors());

        return new Legendre(
                json["polynomial_order"].asInt(),
                nbin
            );
    }
    
    void BasisEvaluator::CoeffInit(void)
    {     
        double coeff_size = 1;
        //Get the size of the number of coefficients + 1
        for(size_t i = 0; i < functions_.size(); ++i)
        {
            coeff_size *= functions_[i]->GetOrder()+1;
        }

        std::vector<int> jdx(functions_.size(), 0);
        Map temp_map(jdx,0.0);
        //Initialize the mapping for the coeff function
        for(size_t i = 0; i < coeff_size; ++i)
        {
            for(size_t j = 0; j < jdx.size(); ++j)
            {
                if(jdx[j] > 0 && jdx[j] % (functions_[j]->GetOrder()+1) == 0)
                {
                    if(j != functions_.size() - 1)
                    { 
                        jdx[j+1]++;
                        jdx[j] = 0;
                    }
                }
                temp_map.map[j] = jdx[j];
                temp_map.value  = 0; 
            }
            coeff_.push_back(temp_map);           
            jdx[0]++;
        }
    }
    
    void BasisEvaluator::BasisInit(void) 
    {
        for (size_t i=0; i<functions_.size(); i++)
        {
            uint poly = functions_[i]->GetOrder()+1;
            uint nbin = functions_[i]->GetBins();
            std::vector<double> dervs(nbin*poly,0);
            std::vector<double> vals(nbin*poly,0);

            for (size_t j=0; j<poly; j++) {
                for(size_t k=0; k<nbin; k++) {
                    double x = ((functions_[i]->GetRange())*k
                            - functions_[i]->GetLower())/nbin + functions_[i]->GetLower();
                    vals[k+j*nbin] = functions_[i]->Evaluate(x,j);
                    dervs[k+j*nbin] = functions_[i]->EvalGrad(x,j);
                }
            }
            BasisLUT TempLUT(vals,dervs);
            lookup_.push_back(TempLUT);
        }
    }
    
    void BasisEvaluator::UpdateBias(Grid<double> *bias, Grid<std::vector<double>> *grad)
    { 
        double basis;
        double temp;
        size_t j = 0;
        int nbins;

        for(Grid<double>::iterator it = bias->begin(); it != bias->end(); ++it, ++j)
        {
            std::vector<double> tmp_grad (functions_.size(),0);
            double tmp_bias = 0;
            for (size_t i = 1; i < coeff_.size(); ++i)
            {
                basis = 1.0;
                for (size_t l = 0; l < functions_.size(); l++)
                {
                    temp = 1.0;
                    //For the gradients we only evaluate the gradient for diagonal terms
                    //Off diagonals are the basis set value
                    for (size_t k = 0; k < functions_.size(); k++)
                    {
                        nbins = bias->GetNumPoints(k);
                        temp *= l == k ?  lookup_[k].derivs[it.index(k) + coeff_[i].map[k]*nbins] * functions_[l]->GetRange() / (bias->GetUpper(l) - bias->GetLower(l))
                                       :  lookup_[k].values[it.index(k) + coeff_[i].map[k]*nbins];
                    }
                    tmp_grad[l] -= coeff_[i].value * temp;
                    nbins = bias->GetNumPoints(l);
                    basis *= lookup_[l].values[it.index(l) + coeff_[i].map[l]*nbins];
                }
                //Update the bias values
                tmp_bias += coeff_[i].value * basis;
                //Store the gradient values
            }
            grad->at(it.coordinates()) = tmp_grad;
            *it = tmp_bias;
        }
    }

    //Calculates the inner product of the basis set and the biased histogram
    //This function then returns the coefficients from this calculation
    double BasisEvaluator::UpdateCoeff(const std::vector<double> &array, Grid<uint> *hist)
    {
        double coeffTemp = 0;
        double sum = 0;

        for(auto& coeff : coeff_)
        {
            // The method uses a standard integration
            size_t j = 0;
            coeffTemp = coeff.value;
            coeff.value = 0.0;
            for(Grid<uint>::iterator it2 = hist->begin(); it2 != hist->end(); ++it2, ++j)
            { 

                double weight = std::pow(2.0,functions_.size());

                // This adds in a trap-rule type weighting which lowers error significantly at the boundaries
                for(size_t k = 0; k < functions_.size(); ++k)
                {
                    if( it2.index(k) == 0 ||
                        it2.index(k) == hist->GetNumPoints(k)-1)
                        weight /= 2.0;
                }
                /*The numerical integration of the biased histogram across the entirety of CV space
                 *All calculations include the normalization as well
                 */
                double basis = 1.0;

                for(size_t l = 0; l < functions_.size(); l++)
                {
                    int nbins = hist->GetNumPoints(l);
                    basis *= lookup_[l].values[it2.index(l) + coeff.map[l]*nbins]*functions_[l]->GetRange()/nbins;
                    //Normalize the values by the associated value
                    basis *= functions_[l]->GetNorm(coeff.map[l])*functions_[l]->Weight(it2.coordinate(l));
                }
                coeff.value += basis * array[j] * weight/std::pow(2.0,functions_.size());
            }
            coeffTemp -= coeff.value;
            //std::cout<<coeffTemp<<std::endl;
            sum += fabs(coeffTemp);
        }
        return sum;
    }
};
