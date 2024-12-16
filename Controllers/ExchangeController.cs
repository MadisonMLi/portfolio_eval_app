using Microsoft.AspNetCore.Mvc;
using System.Text.Json;

namespace Homework6;

[ApiController]
[Route("[controller]")]
public class Exchangeontroller : ControllerBase
{

    [HttpPost("/Exchange")]
    public ActionResult<Exchange> PostExchange([FromBody] Exchange exchange)
    {
        Console.WriteLine("Posted");
        modifydb.modifyexchange(exchange.Id, exchange.Name, exchange.Symbol);
        return Ok(exchange);
    }

    [HttpGet("/Exchange")]
    public ActionResult<Exchange> GetExchange()
    {
        Console.WriteLine("Getting");
        FinanceContext db = new FinanceContext();
        Console.WriteLine(db);
        return Ok(db.Exchanges.ToArray());
    }
}
