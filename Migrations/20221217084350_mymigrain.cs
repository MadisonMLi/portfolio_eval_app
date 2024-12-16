using System;
using Microsoft.EntityFrameworkCore.Migrations;
using Npgsql.EntityFrameworkCore.PostgreSQL.Metadata;

#nullable disable

namespace Homework6.Migrations
{
    /// <inheritdoc />
    public partial class mymigrain : Migration
    {
        /// <inheritdoc />
        protected override void Up(MigrationBuilder migrationBuilder)
        {
            migrationBuilder.CreateTable(
                name: "Exchange",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false)
                        .Annotation("Npgsql:ValueGenerationStrategy", NpgsqlValueGenerationStrategy.IdentityByDefaultColumn),
                    Name = table.Column<string>(type: "text", nullable: false),
                    Symbol = table.Column<string>(type: "text", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_Exchange", x => x.Id);
                });

            migrationBuilder.CreateTable(
                name: "OptionTradeEvaluation",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false)
                        .Annotation("Npgsql:ValueGenerationStrategy", NpgsqlValueGenerationStrategy.IdentityByDefaultColumn),
                    UnrealizedPnl = table.Column<double>(name: "Unrealized_Pnl", type: "double precision", nullable: false),
                    Delta = table.Column<double>(type: "double precision", nullable: false),
                    Gamma = table.Column<double>(type: "double precision", nullable: false),
                    Vega = table.Column<double>(type: "double precision", nullable: false),
                    Rho = table.Column<double>(type: "double precision", nullable: false),
                    Theta = table.Column<double>(type: "double precision", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_OptionTradeEvaluation", x => x.Id);
                });

            migrationBuilder.CreateTable(
                name: "Trade",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false)
                        .Annotation("Npgsql:ValueGenerationStrategy", NpgsqlValueGenerationStrategy.IdentityByDefaultColumn),
                    quantity = table.Column<double>(type: "double precision", nullable: false),
                    financialinstrumentid = table.Column<int>(type: "integer", nullable: false),
                    tradeprice = table.Column<double>(name: "trade_price", type: "double precision", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_Trade", x => x.Id);
                });

            migrationBuilder.CreateTable(
                name: "Unit",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false)
                        .Annotation("Npgsql:ValueGenerationStrategy", NpgsqlValueGenerationStrategy.IdentityByDefaultColumn),
                    typeUnit = table.Column<string>(type: "text", nullable: false),
                    sizeUnit = table.Column<double>(type: "double precision", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_Unit", x => x.Id);
                });

            migrationBuilder.CreateTable(
                name: "Market",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false)
                        .Annotation("Npgsql:ValueGenerationStrategy", NpgsqlValueGenerationStrategy.IdentityByDefaultColumn),
                    Name = table.Column<string>(type: "text", nullable: false),
                    UnitId = table.Column<int>(type: "integer", nullable: false),
                    ExchangeId = table.Column<int>(type: "integer", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_Market", x => x.Id);
                    table.ForeignKey(
                        name: "FK_Market_Exchange_ExchangeId",
                        column: x => x.ExchangeId,
                        principalTable: "Exchange",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                    table.ForeignKey(
                        name: "FK_Market_Unit_UnitId",
                        column: x => x.UnitId,
                        principalTable: "Unit",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateTable(
                name: "FinancialInstrument",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false)
                        .Annotation("Npgsql:ValueGenerationStrategy", NpgsqlValueGenerationStrategy.IdentityByDefaultColumn),
                    marketid = table.Column<int>(type: "integer", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_FinancialInstrument", x => x.Id);
                    table.ForeignKey(
                        name: "FK_FinancialInstrument_Market_marketid",
                        column: x => x.marketid,
                        principalTable: "Market",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateTable(
                name: "Underlying",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false),
                    Name = table.Column<string>(type: "text", nullable: false),
                    Symbol = table.Column<string>(type: "text", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_Underlying", x => x.Id);
                    table.ForeignKey(
                        name: "FK_Underlying_FinancialInstrument_Id",
                        column: x => x.Id,
                        principalTable: "FinancialInstrument",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateTable(
                name: "Option",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false),
                    underlyingid = table.Column<int>(type: "integer", nullable: false),
                    expirationdate = table.Column<DateTime>(name: "expiration_date", type: "timestamp with time zone", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_Option", x => x.Id);
                    table.ForeignKey(
                        name: "FK_Option_FinancialInstrument_Id",
                        column: x => x.Id,
                        principalTable: "FinancialInstrument",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                    table.ForeignKey(
                        name: "FK_Option_Underlying_underlyingid",
                        column: x => x.underlyingid,
                        principalTable: "Underlying",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateTable(
                name: "Asian",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false),
                    Strike = table.Column<double>(type: "double precision", nullable: false),
                    IsCall = table.Column<bool>(name: "Is_Call", type: "boolean", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_Asian", x => x.Id);
                    table.ForeignKey(
                        name: "FK_Asian_Option_Id",
                        column: x => x.Id,
                        principalTable: "Option",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateTable(
                name: "Barrier",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false),
                    Strike = table.Column<double>(type: "double precision", nullable: false),
                    IsCall = table.Column<bool>(name: "Is_Call", type: "boolean", nullable: false),
                    BarrierLevel = table.Column<double>(name: "Barrier_Level", type: "double precision", nullable: false),
                    KnockType = table.Column<int>(name: "Knock_Type", type: "integer", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_Barrier", x => x.Id);
                    table.ForeignKey(
                        name: "FK_Barrier_Option_Id",
                        column: x => x.Id,
                        principalTable: "Option",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateTable(
                name: "Digital",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false),
                    Strike = table.Column<double>(type: "double precision", nullable: false),
                    IsCall = table.Column<bool>(name: "Is_Call", type: "boolean", nullable: false),
                    Payout = table.Column<double>(type: "double precision", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_Digital", x => x.Id);
                    table.ForeignKey(
                        name: "FK_Digital_Option_Id",
                        column: x => x.Id,
                        principalTable: "Option",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateTable(
                name: "European",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false),
                    Strike = table.Column<double>(type: "double precision", nullable: false),
                    IsCall = table.Column<bool>(name: "Is_Call", type: "boolean", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_European", x => x.Id);
                    table.ForeignKey(
                        name: "FK_European_Option_Id",
                        column: x => x.Id,
                        principalTable: "Option",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateTable(
                name: "Lookback",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false),
                    Strike = table.Column<double>(type: "double precision", nullable: false),
                    IsCall = table.Column<bool>(name: "Is_Call", type: "boolean", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_Lookback", x => x.Id);
                    table.ForeignKey(
                        name: "FK_Lookback_Option_Id",
                        column: x => x.Id,
                        principalTable: "Option",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateTable(
                name: "Range",
                columns: table => new
                {
                    Id = table.Column<int>(type: "integer", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_Range", x => x.Id);
                    table.ForeignKey(
                        name: "FK_Range_Option_Id",
                        column: x => x.Id,
                        principalTable: "Option",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateIndex(
                name: "IX_FinancialInstrument_marketid",
                table: "FinancialInstrument",
                column: "marketid");

            migrationBuilder.CreateIndex(
                name: "IX_Market_ExchangeId",
                table: "Market",
                column: "ExchangeId");

            migrationBuilder.CreateIndex(
                name: "IX_Market_UnitId",
                table: "Market",
                column: "UnitId");

            migrationBuilder.CreateIndex(
                name: "IX_Option_underlyingid",
                table: "Option",
                column: "underlyingid");
        }

        /// <inheritdoc />
        protected override void Down(MigrationBuilder migrationBuilder)
        {
            migrationBuilder.DropTable(
                name: "Asian");

            migrationBuilder.DropTable(
                name: "Barrier");

            migrationBuilder.DropTable(
                name: "Digital");

            migrationBuilder.DropTable(
                name: "European");

            migrationBuilder.DropTable(
                name: "Lookback");

            migrationBuilder.DropTable(
                name: "OptionTradeEvaluation");

            migrationBuilder.DropTable(
                name: "Range");

            migrationBuilder.DropTable(
                name: "Trade");

            migrationBuilder.DropTable(
                name: "Option");

            migrationBuilder.DropTable(
                name: "Underlying");

            migrationBuilder.DropTable(
                name: "FinancialInstrument");

            migrationBuilder.DropTable(
                name: "Market");

            migrationBuilder.DropTable(
                name: "Exchange");

            migrationBuilder.DropTable(
                name: "Unit");
        }
    }
}
