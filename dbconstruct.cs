namespace Homework6;

using System;
using System.Text;
using System.Collections.Generic;
using Microsoft.EntityFrameworkCore;
using Microsoft.EntityFrameworkCore.Design;
using System.ComponentModel.DataAnnotations.Schema; 
using System.Linq;
using Npgsql;

public class modifydb
{
    // migration name is "defaultmigration", dont know if that matters
    //
    public static void modifyexchange(int id, string name, string symbol)
    {
        var db = new FinanceContext();

        db.Add(new Exchange{
                Id = id,
                Name = name,
                Symbol = symbol,
        });
        db.SaveChanges();
    }
    // public static void modifymarkets(int id, string name, Unit our_unit, Exchange our_exhange)
    public static void modifymarkets(int id, string name, int unitid, Unit myunit, int exchangeid, Exchange myexchange)
    {
        var db = new FinanceContext();

        db.Units.Attach(myunit);
        db.Exchanges.Attach(myexchange);

        db.Add(new Market{
                Id = id,
                Name = name,
                UnitId = unitid,
                MyUnit = myunit,
                ExchangeId = exchangeid,
                MyExchange = myexchange
        });
        db.SaveChanges();
    }
    public static void modifyunits(int id, string type_unit, double size_unit)
    {
        var db = new FinanceContext();
        db.Add(new Unit{
                Id = id,
                typeUnit = type_unit,
                sizeUnit = size_unit,
        });
        db.SaveChanges();
    }
    public static void modifyfinancialinstrument(int id, Market market)
    {
        var db = new FinanceContext();
        db.Markets.Attach(market);
        db.Add(new FinancialInstrument{
                Id = id,
                marketid = market.Id,
                MyMarket = market,
        });
        db.SaveChanges();
    }
    // public static void modifyratepoints(int id, DateTime our_date, double our_term, double our_rate, Ratecurve our_ratecurve)
    // {
    //     var db = new FinanceContext();
    //     db.Add(new RatePoint{
    //             Id = id,
    //             measure_date = our_date,
    //             term = our_term,
    //             rate = our_rate,
    //             ratecurve_id = our_ratecurve.Id,
    //     });
    //     db.SaveChanges();
    // }
    // public static void modifyratecurves(int id, string name, string symbol)
    // {
    //     var db = new FinanceContext();
    //     db.Add(new Ratecurve{
    //             Id = id,
    //             name = name,
    //             symbol = symbol,
    //     });
    //     db.SaveChanges();
    // }
    // public static void modifyhistoricalprices(int id,DateTime our_date, string our_type, int types_id)
    // {
    //     var db = new FinanceContext();
    //     db.Add(new Historical_price{
    //             Id = id,
    //             dateofobservation = our_date,
    //             instrument = our_type,
    //             instrument_id= types_id,
    //     });
    //     db.SaveChanges();
    // }
    public static void modifyunderlyings(int id, int marketid, string name, string symbol, Market market)
    {
        var db = new FinanceContext();

        db.Markets.Attach(market);

        db.Add(new Underlying{
                Id = id,
                Name = name,
                Symbol = symbol,
                marketid = marketid,
                MyMarket = market,
        });
        db.SaveChanges();
    }
    public static void modifyoptions(Underlying our_underlying, DateTime our_date)
    { 
        var db = new FinanceContext();
        db.Add(new Option{
                Id = our_underlying.Id,
                expiration_date = our_date,
        });
        db.SaveChanges();
    }
    public static void modifyeuropean(Market market, Underlying underlying, DateTime expiration_date, double strike, Boolean iscall)
    { 
        var db = new FinanceContext();
        db.Markets.Attach(market);
        db.Underlyings.Attach(underlying);

        db.Add(new European{
            marketid = market.Id,
            MyMarket = market,
            underlyingid = underlying.Id,
            expiration_date = expiration_date,
            Strike = strike,
            Is_Call = iscall,
        });
        db.SaveChanges();
    }
    public static void modifydigital(Market market, Underlying underlying, DateTime expiration_date, double strike, Boolean iscall, double payout)
    {
        var db = new FinanceContext();
        db.Markets.Attach(market);
        db.Underlyings.Attach(underlying);
        db.Add(new Digital{
                Strike = strike,
                Is_Call = iscall,
                marketid = market.Id,
                MyMarket = market,
                underlyingid = underlying.Id,
                expiration_date = expiration_date,
                Payout = payout
        });
        db.SaveChanges();
    }
    public static void modifyasian(Market market, Underlying underlying, DateTime expiration_date, double strike, Boolean iscall)
    { 
        var db = new FinanceContext();
        db.Markets.Attach(market);
        db.Underlyings.Attach(underlying);

        db.Add(new Asian{
                Strike = strike,
                Is_Call = iscall,
                marketid = market.Id,
                MyMarket = market,
                underlyingid = underlying.Id,
                expiration_date = expiration_date,
        });
        db.SaveChanges();
    }
    public static void modifyrange(Market market, Underlying underlying, DateTime expiration_date)
    { 
        var db = new FinanceContext();
        db.Markets.Attach(market);
        db.Underlyings.Attach(underlying);

        db.Add(new Range{
                marketid = market.Id,
                MyMarket = market,
                underlyingid = underlying.Id,
                expiration_date = expiration_date,
        });
        db.SaveChanges();
    }
    public static void modifylookback(Market market, Underlying underlying, DateTime expiration_date, double strike, Boolean iscall)
    { 
        var db = new FinanceContext();
        db.Markets.Attach(market);
        db.Underlyings.Attach(underlying);
        db.Add(new Lookback{
                Strike = strike,
                Is_Call = iscall,
                marketid = market.Id,
                MyMarket = market,
                underlyingid = underlying.Id,
                expiration_date = expiration_date,
        });
        db.SaveChanges();
    }

    //dont need to modify range, it has no real properties
    public static void modifybarrier(Market market, Underlying underlying, DateTime expiration_date, double strike, Boolean iscall, double our_barrier_level, int knock_num)
    {
        var db = new FinanceContext();
        db.Markets.Attach(market);
        db.Underlyings.Attach(underlying);
        db.Add(new Barrier{
                Strike = strike,
                Is_Call = iscall,
                marketid = market.Id,
                MyMarket = market,
                underlyingid = underlying.Id,
                expiration_date = expiration_date,
                Barrier_Level = our_barrier_level,
                Knock_Type = knock_num,
        });
        db.SaveChanges();
    }
    public static void modifytrades(int id, double myquantity, int myfinancialinstrumentid, double mytrade_price)
    { 
        /*
            public int Id {get;set;}
            public double quantity {get;set;} // change to quantity, dont need is buy
            public int underlyingid {get;set;}
            public string instrument_type {get;set;} // option, future, underlying, etc.
            public double Trade_Price {get;set;} // price it was traded at
        */
        var db = new FinanceContext();
        db.Add(new Trade{
                Id = id,
                quantity = myquantity,
                financialinstrumentid = myfinancialinstrumentid,
                trade_price = mytrade_price,
        });
        db.SaveChanges();
    }
    public static void modifyoptiontradeevaluations(int id,double our_pnl, double our_delta, double our_gamma, double our_vega, double our_rho, double our_theta)
    { 
        /*
        public int Id {get;set;}
    
        // double standard_error {get; set;}
        // Monte Carlo evals. an instrument; not a trade.
        public double Unrealized_Pnl {get; set;} // different from market value; refer to trade price
        public double Delta {get; set;}
        public double Gamma {get; set;}
        public double Vega {get; set;}
        public double Rho {get; set;}
        public double Theta {get; set;}
        */
        var db = new FinanceContext();
        db.Add(new Option_Trade_Evaluation{
               Id = id,
               Unrealized_Pnl = our_pnl,
               Delta = our_delta,
               Gamma = our_gamma,
               Vega = our_vega,
               Rho = our_rho,
               Theta = our_theta,
        });
        db.SaveChanges();
    }
}
public class FinanceContext : DbContext{
    public DbSet<Exchange> Exchanges { get; set;}
    public DbSet<Market> Markets {get;set;}
    public DbSet<Unit> Units {get;set;}
    public DbSet<FinancialInstrument> FinancialInstruments {get;set;}
    // public DbSet<RatePoint> RatePoints {get;set;}
    // public DbSet<Ratecurve> RateCurves {get;set;}
    // public DbSet<Historical_price> HistoricalPrices {get;set;}
    public DbSet<Underlying> Underlyings {get; set;}
    public DbSet<Option> Options {get;set;}
    public DbSet<European> EuropeanOptions {get;set;}
    public DbSet<Digital> DigitalOptions {get;set;}
    public DbSet<Asian> AsianOptions {get;set;}
    public DbSet<Lookback> LookbackOptions {get;set;}
    public DbSet<Range> RangeOptions {get;set;}
    public DbSet<Barrier> BarrierOptions{get;set;}
    public DbSet<Trade> Trades{get;set;}
    public DbSet<Option_Trade_Evaluation> OptionTradeEvaluations{get;set;}

    protected override void OnConfiguring(DbContextOptionsBuilder optionsBuilder)
                    => optionsBuilder.UseNpgsql("host=localhost; Database=postgres; Username=postgres; Password=secret");
}

[Table("Exchange")]
public class Exchange
{
    public int Id {get;set;}
    public string Name {get;set;}
    public string Symbol {get;set;}
    //public List<Market> markets {get; set;}
}

[Table("Market")]
public class Market
{
    public int Id{get;set;}
    public string Name {get;set;}
    public int UnitId {get;set;}
    public Unit ? MyUnit {get;set;}
    public int ExchangeId {get;set;}
    public Exchange ? MyExchange {get;set;}
}

[Table("Unit")]
public class Unit
{   public int Id{get;set;}
    public string typeUnit {get;set;}
    public double sizeUnit {get;set;}
}

// [Table("RatePoint")]
// public class RatePoint
// {
//     public int Id{get;set;}
//     public DateTime measure_date {get;set;} // put into historical price?
//     public double term {get;set;}
//     public double rate {get;set;}
//     public int ratecurve_id {get;set;} //this identifies which rate curve it belongs to
// }

// [Table("RateCurve")]
// public class Ratecurve
// {
//     public int Id{get;set;}
//     public string name {get;set;}
//     public string symbol {get;set;}
//     public List<RatePoint> ? ratepoints {get;set;}
// }

// [Table("HistoricalPrice")] //I honestly dont think we'll be using this much, I cant truly think of a reason to
// public class Historical_price
// {
//     public int Id{get;set;}
//     public List<double> prices {get; set;} // this means each of our prices need to be linked... does this need to exist?
//     public DateTime dateofobservation {get; set;}
//     public string instrument {get; set;} // this is either underlying, option, or future
//     public int instrument_id {get; set;} // price_of to go to the proper table and price_of_id to go to proper option, underlying, or future
//     // id, refers to inst, date of obser, price(s)
// }

[Table("FinancialInstrument")]
public class FinancialInstrument
{   
    public int Id{get;set;}
    public int marketid {get;set;}
    public Market ? MyMarket {get;set;}
}

[Table("Underlying")]
public class Underlying : FinancialInstrument
{
    // market underlying_market{get;set;}
    // public int Id{get;set;}
    // public int marketid {get; set;}
    public string Name {get;set;}
    public string Symbol {get;set;}
    //historical_close_prices prices {get;set;}
}

[Table("Option")]
public class Option: FinancialInstrument
{
    public int underlyingid {get;set;}
    public Underlying ? underlying {get;set;}
    public DateTime expiration_date {get; set;}
}

[Table("European")]
public class European : Option
{
    public double Strike {get; set;}
    public bool Is_Call {get; set;}
}

[Table("Digital")]

public class Digital : Option
{
    public double Strike {get; set;}
    public Boolean Is_Call {get; set;}
    public double Payout {get; set;}
}

[Table("Asian")]
public class Asian : Option
{
    public double Strike {get; set;}
    public Boolean Is_Call {get; set;}
}
[Table("Range")]
public class Range : Option
{
    // Nothing here, dawg.
}

[Table("Lookback")]
public class Lookback : Option
{
    public double Strike {get; set;}
    public Boolean Is_Call {get; set;}
}

[Table("Barrier")]
public class Barrier : Option
{
    public double Strike {get; set;}
    public Boolean Is_Call {get; set;}
    public double Barrier_Level {get; set;}
    public int Knock_Type {get; set;}
}

[Table("Trade")]
public class Trade
{  
    public int Id {get;set;}
    public double quantity {get;set;} // change to quantity, dont need is buy
    public int financialinstrumentid {get;set;}
    // public string instrument_type {get;set;} // option, future, underlying, etc.
    public double trade_price {get;set;} // price it was traded at
}

[Table("OptionTradeEvaluation")]
public class Option_Trade_Evaluation
{
    public int Id {get;set;}
    
    // double standard_error {get; set;}
    // Monte Carlo evals. an instrument; not a trade.
    public double Unrealized_Pnl {get; set;} // different from market value; refer to trade price
    public double Delta {get; set;}
    public double Gamma {get; set;}
    public double Vega {get; set;}
    public double Rho {get; set;}
    public double Theta {get; set;}
}
