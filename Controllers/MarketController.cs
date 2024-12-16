using Microsoft.AspNetCore.Mvc;
using System.Text.Json;

namespace Homework6;

[ApiController]
[Route("[controller]")]
public class MarketController : ControllerBase
{

    [HttpPost("/Market")]
    public ActionResult<Market> PostMarket([FromBody] Market market)
    {
        Console.WriteLine("API Post");
        FinanceContext db = new FinanceContext();

        Unit myunit = db.Units.Single(x => x.Id == market.UnitId);
        Exchange myexchange = db.Exchanges.Single(x => x.Id == market.ExchangeId);

        modifydb.modifymarkets(market.Id, market.Name, market.UnitId, myunit, market.ExchangeId, myexchange);

        return Ok(market);
    }

    [HttpGet("/Market")]
    public ActionResult<Market> GetMarket()
    {
        Console.WriteLine("API Get");
        FinanceContext db = new FinanceContext();
        Console.WriteLine(db);
        return Ok(db.Markets.ToArray());
    }
}
