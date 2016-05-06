# ----------------------------------------------------------------------------
#
#  carray/math/gsl.rb
#
#  This file is part of Ruby/CArray extension library.
#  You can redistribute it and/or modify it under the terms of
#  the GNU General Public License (GPL).
#
#  Copyright (C) 2005-2008  Hiroki Motoyoshi
#
# ----------------------------------------------------------------------------

require 'gsl'
require 'carray'
require 'carray/carray_gsl'

module CA::GSL
end

[ 
 GSL::Vector,      GSL::Vector::Int,      GSL::Vector::Complex,
 GSL::Vector::Col, GSL::Vector::Int::Col, GSL::Vector::Complex::Col,
 GSL::Matrix,      GSL::Matrix::Int,      GSL::Matrix::Complex,
].each do |klass|
  klass.class_eval {
    def to_ca
      return self.ca.to_ca
    end
  }
end

class CArray

  [
    ["airy_Ai", []],
    ["airy_Bi", []],
    ["airy_Ai_scaled", []],
    ["airy_Bi_scaled", []],
    ["airy_deriv_Ai", []],
    ["airy_deriv_Bi", []],
    ["airy_deriv_Ai_scaled", []],
    ["airy_deriv_Bi_scaled", []],
  ].each do |name,parms|
    parms1 = parms + ["mode=GSL::PREC_DOUBLE"]
    parms2 = parms + ["mode"]
    eval %{
      def #{name} (#{parms1.join(",")})
        return self.#{name}_(#{parms2.join(",")})
      end
    }
  end

  [
   ["random_gaussian", ["sigma"]],
   ["random_ugaussian", []],
   ["random_gaussian_tail", ["a", "sigma"]],
   ["random_ugaussian_tail", ["a"]],
   ["random_exponential", ["mu"]],
   ["random_laplace", ["a"]],
   ["random_exppaw", ["a", "b"]],
   ["random_cauchy", ["a"]],
   ["random_rayleigh", ["sigma"]],
   ["random_rayleigh_tail", ["a", "sigma"]],
   ["random_landau", []],
   ["random_levy", ["c", "alpha"]],
   ["random_levy_skew", ["c", "alpha", "beta"]],
   ["random_gamma", ["a", "b"]],
   ["random_erlang", ["a", "n"]],
   ["random_flat", ["a", "b"]],
   ["random_lognormal", ["zeta", "sigma"]],
   ["random_chisq", ["nu"]],
   ["random_fdist", ["nu1", "nu2"]],
   ["random_tdist", ["nu"]],
   ["random_beta", ["a", "b"]],
   ["random_logistic", ["a"]],
   ["random_pareto", ["a", "b"]],
   ["random_weibull", ["a", "b"]],
   ["random_gumbel1", ["a", "b"]],
   ["random_gumbel2", ["a", "b"]],
   ["random_poisson", ["mu"]],
   ["random_bernoulli", ["p"]],
   ["random_binomial", ["p", "n"]],
   ["random_negative_binomial", ["p", "n"]],
   ["random_pascal", ["p", "n"]],
   ["random_geometric", ["p"]],
   ["random_hypergeometric", ["n1", "n2", "t"]],
   ["random_logarithmic", ["p"]],
  ].each do |name, parms|
    eval %{
      def #{name} (#{parms.join(",")})
        return CAMath.#{name}(#{(["self"]+parms).join(",")}).to_type(data_type)
      end
      def #{name}! (#{parms.join(",")})
        return self.asign{ CAMath.#{name}(#{(["self"]+parms).join(",")}) }
      end
    }
  end

  def lu_decomp
    attach {
      lu, perm, sign = *GSL::Linalg::LU.decomp(self.gm)
      out = lu.ca.m
      out.attribute[:A]    = self.gm
      out.attribute[:LU]   = lu
      out.attribute[:perm] = perm
      out.attribute[:sign] = sign
      return out
    }
  end

  def lu_solve (b)
    if self.attribute[:LU]
      CArray.attach(self, b) {
        lu   = self.attribute[:LU]
        perm = self.attribute[:perm]
        return GSL::Linalg::LU.solve(lu, perm, b.gv).ca.v
      }
    else
      CArray.attach(self, b) {
        return GSL::Linalg::LU.solve(self.gm, b.gv).ca.v
      }
    end
  end

  def lu_refine (b, x)
    CArray.attach(self, b, x) {
      a    = self.attribute[:A]
      lu   = self.attribute[:LU]
      perm = self.attribute[:perm]
      return GSL::Linalg::LU.refine(a, lu, perm, b.gv, x.gv).first.ca.v
    }
  end

  def lu_invert
    attach {
      return GSL::Linalg::LU.invert(self.gm).ca.m
    }
  end

  def lu_det
    attach {
      return GSL::Linalg::LU.det(self.gm)
    }
  end

  def qr_decomp
    attach {
      qr, tau = *GSL::Linalg::QR.decomp(self.gm)
      out = qr.ca.m
      out.attribute[:qr] = qr
      out.attribute[:tau] = tau
      return out
    }
  end

  def qr_solve (b)
    CArray.attach(self, b) {
      return GSL::Linalg::QR.solve(self.gm, b.gv).ca.v
    }
  end

  def qr_unpack 
    attach {
      qr  = self.attribute[:qr]
      tau = self.attribute[:tau]
      q, r = *GSL::Linalg::QR.unpack(qr, tau)
      return q.ca.m, r.ca.m
    }
  end

  def sv_decomp
    attach {
      u, v, s = *GSL::Linalg::SV.decomp(self.gm)
      return u.ca.m, v.ca.m, s.ca.v
    }
  end

  def sv_solve (b)
    CArray.attach(self, b) {
      return GSL::Linalg::SV.solve(self.gm, b.gv).ca.v
    }
  end

  def cholesky_decomp
    attach {
      c = GSL::Linalg::Cholesky.decomp(self.gm)
      out = c.ca.m
      out.attribute[:cholesky] = c
      return out
    }
  end

  def cholesky_solve (b)
    CArray.attach(self, b) {
      return GSL::Linalg::Cholesky.solve(self.gm, b.gv).ca.v
    }
  end

  def hh_solve (b)
    CArray.attach(self, b) {
      return GSL::Linalg::HH.solve(self.gm, b.gv).ca.v
    }
  end

  def eigen_symm 
    attach {
      return GSL::Eigen.symm(self.gm).ca.v
    }
  end
  
  def eigen_symmv
    attach {
      e, v = *GSL::Eigen.symmv(self.gm)
      return e.ca.v, v.ca.m
    }
  end
  
  def eigen_herm
    attach {
      return GSL::Eigen.herm(self.gm).ca.v
    }
  end
  
  def eigen_hermv
    attach {
      e, v = *GSL::Eigen.hermv(self.gm)
      return e.ca.v, v.ca.m
    }
  end
  
  def eigen_nonsymm 
    attach {
      return GSL::Eigen.nonsymm(self.gm).ca.v
    }
  end
  
  def eigen_nonsymm_Z 
    attach {
      e, z = *GSL::Eigen.nonsymm_Z(self.gm)
      return e.ca.v, z.ca.m
    }
  end
  
  def eigen_nonsymmv 
    attach {
      e, v = *GSL::Eigen.nonsymmv(self.gm)
      return e.ca.v, v.ca.m
    }
  end

  def fit_linear (x)
    CArray.attach(self, x) {
      c0, c1, c00, c01, c11, chi2, status = 
                         *GSL::Fit.linear(x.gv, self.gv)
      dof   = self.elements - 2
      covar = CA_DOUBLE([[c00, c01],[c01, c11]])
      err   = Array.new(2){|i| Math::sqrt(chi2/dof*covar[i,i]) }
      return [c0, c1], err, chi2, dof, covar
    }
  end

  def fit_wlinear (x, w)
    CArray.attach(self, x, w) {
      c0, c1, c00, c01, c11, chi2, status = 
                         *GSL::Fit.wlinear(x.gv, w.gv, self.gv)
      dof   = self.elements - 2
      covar = CA_DOUBLE([[c00, c01],[c01, c11]])
      err   = Array.new(2){|i| Math::sqrt(chi2/dof*covar[i,i]) }
      return [c0, c1], err, chi2, dof, covar
    }
  end

  def fit_nonlinear (xx, sigma, procf, procdf, parms, errs)
    n  = self.elements
    np = parms.size
    parms = GSL::Vector[*parms]
    gf = Proc.new { |a, x, y, sigma, f|
      f.ca[] = (procf.call(x, a.to_a) - y) / sigma
    }
    gdf = Proc.new { |a, x, y, sigma, jac| 
      cjac = jac.ca
      df   = procdf.call(x, a.to_a)
      a.size.times do |i|
        cjac[nil, i] = df[i] / sigma
      end
    }
    func   = GSL::MultiFit::Function_fdf.alloc(gf, gdf, np)
    func.set_data(xx, self, sigma)
    lmsder = GSL::MultiFit::FdfSolver::LMSDER
    solver = GSL::MultiFit::FdfSolver.alloc(lmsder, n, np)
    solver.set(func, parms)

    iter = 0
    begin
      iter += 1
      status = solver.iterate
      status = solver.test_delta(*errs)
    end while status == GSL::CONTINUE and iter < 500

    coef  = solver.position.to_a
    chi2  = GSL::pow_2(solver.f.dnrm2)
    dof   = n - np
    covar = solver.covar(0.0).ca
    err   = Array.new(np){|i| Math::sqrt(chi2/dof*covar[i,i]) }
    return coef, err, chi2, dof, covar
  end

  def fit_builtin (*args)
    if args[1].is_a?(String)
      x, func, guess = *args
      CArray.attach(self, x) {
        if guess
          coef, err, chi2, dof =
            *GSL::MultiFit::FdfSolver.fit(x.gv, self.gv, func, guess)
        else
          coef, err, chi2, dof =
            *GSL::MultiFit::FdfSolver.fit(x.gv, self.gv, func)
        end
        return coef.to_a, err.to_a, chi2, dof
      }
    elsif args[2].is_a?(String)
      x, w, func, guess = *args
      CArray.attach(self, x, w) {
        if guess
          coef, err, chi2, dof =
            *GSL::MultiFit::FdfSolver.fit(x.gv, w.gv, self.gv, func, guess)
        else
          coef, err, chi2, dof =
            *GSL::MultiFit::FdfSolver.fit(x.gv, w.gv, self.gv, func)
        end
        return coef.to_a, err.to_a, chi2, dof
      }
    else
      raise ArgumentError, "invalid built-in function specifier"
    end
  end

  def stats_mean 
    attach {
      return GSL::Stats.mean(self.gv)
    }
  end

  def stats_variance 
    attach {
      return GSL::Stats.variance(self.gv)
    }
  end

  def stats_variance_m (mean)
    attach {
      return GSL::Stats.variance_m(self.gv, mean)
    }
  end

  def stats_sd 
    attach {
      return GSL::Stats.sd(self.gv)
    }
  end

  def stats_sd_m (mean)
    attach {
      return GSL::Stats.sd_m(self.gv, mean)
    }
  end

  def stats_tss 
    attach {
      return GSL::Stats.tss(self.gv)
    }
  end

  def stats_tss_m (mean)
    attach {
      return GSL::Stats.tss_m(self.gv, mean)
    }
  end

  def stats_variance_with_fixed_mean (mean)
    attach {
      return GSL::Stats.variance_with_fixed_mean(self.gv, mean)
    }
  end

  def stats_sd_with_fixed_mean (mean)
    attach {
      return GSL::Stats.sd_with_fixed_mean(self.gv, mean)
    }
  end

  def stats_absdev 
    attach {
      return GSL::Stats.absdev(self.gv)
    }
  end

  def stats_absdev_m (mean)
    attach {
      return GSL::Stats.absdev_m(self.gv, mean)
    }
  end

  def stats_skew 
    attach {
      return GSL::Stats.skew(self.gv)
    }
  end

  def stats_skew_m_sd (mean, sd)
    attach {
      return GSL::Stats.skew(self.gv, mean, sd)
    }
  end

  def stats_kurtosis 
    attach {
      return GSL::Stats.kurtosis(self.gv)
    }
  end

  def stats_kurtosis_m_sd (mean, sd)
    attach {
      return GSL::Stats.kurtosis(self.gv, mean, sd)
    }
  end

  def stats_lag1_autocorrelation 
    attach {
      return GSL::Stats.lag1_autocorrelation(self.gv)
    }
  end

  def stats_lag1_autocorrelation_m (mean)
    attach {
      return GSL::Stats.lag1_autocorrelation(self.gv, mean)
    }
  end


  def stats_wmean (w)
    CArray.attach(self, w) {
      return GSL::Stats.wmean(self.gv, w.gv)
    }
  end

  def stats_wvariance (w)
    CArray.attach(self, w) {
      return GSL::Stats.wvariance(self.gv, w.gv)
    }
  end

  def stats_wvariance_m (w, mean)
    CArray.attach(self, w) {
      return GSL::Stats.wvariance_m(self.gv, w.gv, mean)
    }
  end

  def stats_wsd (w)
    CArray.attach(self, w) {
      return GSL::Stats.wsd(self.gv, w.gv)
    }
  end

  def stats_wsd_m (w, mean)
    CArray.attach(self, w) {
      return GSL::Stats.wsd_m(self.gv, w.gv, mean)
    }
  end

  def stats_wtss (w)
    CArray.attach(self, w) {
      return GSL::Stats.wtss(self.gv, w.gv)
    }
  end

  def stats_wtss_m (w, mean)
    CArray.attach(self, w) {
      return GSL::Stats.wtss_m(self.gv, w.gv, mean)
    }
  end

  def stats_wvariance_with_fixed_mean (w, mean)
    CArray.attach(self, w) {
      return GSL::Stats.wvariance_with_fixed_mean(self.gv, w.gv, mean)
    }
  end

  def stats_wsd_with_fixed_mean (w, mean)
    CArray.attach(self, w) {
      return GSL::Stats.wsd_with_fixed_mean(self.gv, w.gv, mean)
    }
  end

  def stats_wabsdev (w)
    CArray.attach(self, w) {
      return GSL::Stats.wabsdev(self.gv, w.gv)
    }
  end

  def stats_wabsdev_m (w, mean)
    CArray.attach(self, w) {
      return GSL::Stats.wabsdev_m(self.gv, w.gv, mean)
    }
  end

  def stats_wskew (w)
    CArray.attach(self, w) {
      return GSL::Stats.wskew(self.gv, w.gv)
    }
  end

  def stats_wskew_m_sd (w, mean, sd)
    CArray.attach(self, w) {
      return GSL::Stats.wskew(self.gv, w.gv, mean, sd)
    }
  end

  def stats_wkurtosis (w)
    CArray.attach(self, w) {
      return GSL::Stats.wkurtosis(self.gv, w.gv)
    }
  end

  def stats_wkurtosis_m_sd (w, mean, sd)
    CArray.attach(self, w) {
      return GSL::Stats.wkurtosis(self.gv, w.gv, mean, sd)
    }
  end

  def stats_max 
    attach {
      return GSL::Stats.max(self.gv)
    }
  end

  def stats_min 
    attach {
      return GSL::Stats.min(self.gv)
    }
  end

  def stats_minmax
    attach {
      return GSL::Stats.minmax(self.gv)
    }
  end

  def stats_max_index
    attach {
      return GSL::Stats.max_index(self.gv)
    }
  end

  def stats_min_index
    attach {
      return GSL::Stats.min_index(self.gv)
    }
  end

  def stats_minmax_index
    attach {
      return GSL::Stats.minmax_index(self.gv)
    }
  end
  
  def stats_median_from_sorted_data
    attach {
      return GSL::Stats.median_from_sorted_data(self.gv)
    }
  end

  def stats_quantile_from_sorted_data
    attach {
      return GSL::Stats.quantile_from_sorted_data(self.gv)
    }
  end

end

module CAMath

  include GSL::CONST::MKSA
  include GSL::CONST::NUM

  module_function

  def coupling_9j (*argv)
    return GSL.Sf.coupling_9j(*argv)
  end

  def lu_solve (lu, b)
    CArray.attach(lu, b) {
      perm = lu.attribute[:perm]
      return GSL::Linalg::LU.solve(lu.gm, perm, b.gv).ca.v
    }
  end

  def qr_solve (qr, b)
    CArray.attach(qr, b) {
      qrm = qr.attribute[:qr]
      tau = qr.attribute[:tau]
      return GSL::Linalg::QR.solve(qrm, tau, b.gv).ca.v
    }
  end

  def sv_solve (u, v, s, b)
    CArray.attach(u, v, s, b) {
      return GSL::Linalg::SV.solve(u.gm, v.gm, s.gv, b.gv).ca.v
    }
  end

  def cholesky_solve (c, b)
    CArray.attach(c, b) {
      cholesky = c.attribute[:cholesky]
      return GSL::Linalg::Cholesky.solve(cholesky, b.gv).ca.v
    }
  end

  def stats_covariance (v1, v2)
    CArray.attach(v1, v2) {
      return GSL::Stats.covariance(v1.gv, v2.gv)
    }
  end

  def stats_covariance_m (v1, v2, mean1, mean2)
    CArray.attach(v1, v2) {
      return GSL::Stats.covariance_m(v1.gv, v2.gv, mean1, mean2)
    }
  end

  def stats_correlation (v1, v2)
    CArray.attach(v1, v2) {
      return GSL::Stats.correlation(v1.gv, v2.gv)
    }
  end

  def graph (*args)
    if args.last.is_a?(String)
      options = args.pop
    else
      options = nil
    end
    args.map! do |arg|
      case arg
      when CArray
        arg.to_gv
      when Array
        arg.map{|a| a.is_a?(CArray) ? a.to_gv : a }
      else
        arg
      end
    end
    return GSL.graph(*args)
  end

end

class CAMatrix 

  def * (other)
    case other
    when CAVector
      CArray.attach(self, other) {
        return (self.gm * other.gv).ca.v      
      }
    when CAMatrix
      CArray.attach(self, other) {
        return (self.gm * other.gm).ca.m      
      }
    when CArray
      CArray.attach(self, other) {
        return (self[] * other).m
      }
    else
      attach {
        return (self.gm * other).ca.m      
      }
    end
  end

  def lower
    attach {
      return self.gm.lower.ca.m
    }
  end

  def upper
    attach {
      return self.gm.lower.ca.m
    }
  end

  def diagonal
    attach {
      return self.gm.diagonal.ca.v
    }
  end

end

class CAVector

  def * (other)
    case other
    when CAVector
      CArray.attach(self, other) {
        out = self.gv * other.gv
        case out
        when GSL::Matrix
          return out.ca.m
        when GSL::Vector
          return out.ca.v
        else
          return out
        end
      }
    when CAMatrix
      CArray.attach(self, other) {
        return (self.gv * other.gm).ca.v
      }
    when CArray
      CArray.attach(self, other) {
        return (self[] * other).v
      }
    else
      attach {
        return (self.gv * other).ca.v   
      }
    end
  end

  alias t_orig t

  def t
    case rank
    when 1
      return self.reshape(elements,1).v
    else
      return t_orig.v
    end
  end

  def norm
    attach {
      return self.gv.nrm2
    }
  end

  def normalize (nrm=1.0)
    attach {
      return self.gv.normalize(nrm).ca.v
    }
  end

  def normalize!
    attach {
      self.gv.normalize!(nrm)
    }
    return self
  end

end

class CArray
  def vt
    return v.t
  end
end
