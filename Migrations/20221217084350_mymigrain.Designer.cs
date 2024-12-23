﻿// <auto-generated />
using System;
using Homework6;
using Microsoft.EntityFrameworkCore;
using Microsoft.EntityFrameworkCore.Infrastructure;
using Microsoft.EntityFrameworkCore.Migrations;
using Microsoft.EntityFrameworkCore.Storage.ValueConversion;
using Npgsql.EntityFrameworkCore.PostgreSQL.Metadata;

#nullable disable

namespace Homework6.Migrations
{
    [DbContext(typeof(FinanceContext))]
    [Migration("20221217084350_mymigrain")]
    partial class mymigrain
    {
        /// <inheritdoc />
        protected override void BuildTargetModel(ModelBuilder modelBuilder)
        {
#pragma warning disable 612, 618
            modelBuilder
                .HasAnnotation("ProductVersion", "7.0.0")
                .HasAnnotation("Relational:MaxIdentifierLength", 63);

            NpgsqlModelBuilderExtensions.UseIdentityByDefaultColumns(modelBuilder);

            modelBuilder.Entity("Homework6.Exchange", b =>
                {
                    b.Property<int>("Id")
                        .ValueGeneratedOnAdd()
                        .HasColumnType("integer");

                    NpgsqlPropertyBuilderExtensions.UseIdentityByDefaultColumn(b.Property<int>("Id"));

                    b.Property<string>("Name")
                        .IsRequired()
                        .HasColumnType("text");

                    b.Property<string>("Symbol")
                        .IsRequired()
                        .HasColumnType("text");

                    b.HasKey("Id");

                    b.ToTable("Exchange");
                });

            modelBuilder.Entity("Homework6.FinancialInstrument", b =>
                {
                    b.Property<int>("Id")
                        .ValueGeneratedOnAdd()
                        .HasColumnType("integer");

                    NpgsqlPropertyBuilderExtensions.UseIdentityByDefaultColumn(b.Property<int>("Id"));

                    b.Property<int>("marketid")
                        .HasColumnType("integer");

                    b.HasKey("Id");

                    b.HasIndex("marketid");

                    b.ToTable("FinancialInstrument");

                    b.UseTptMappingStrategy();
                });

            modelBuilder.Entity("Homework6.Market", b =>
                {
                    b.Property<int>("Id")
                        .ValueGeneratedOnAdd()
                        .HasColumnType("integer");

                    NpgsqlPropertyBuilderExtensions.UseIdentityByDefaultColumn(b.Property<int>("Id"));

                    b.Property<int>("ExchangeId")
                        .HasColumnType("integer");

                    b.Property<string>("Name")
                        .IsRequired()
                        .HasColumnType("text");

                    b.Property<int>("UnitId")
                        .HasColumnType("integer");

                    b.HasKey("Id");

                    b.HasIndex("ExchangeId");

                    b.HasIndex("UnitId");

                    b.ToTable("Market");
                });

            modelBuilder.Entity("Homework6.Option_Trade_Evaluation", b =>
                {
                    b.Property<int>("Id")
                        .ValueGeneratedOnAdd()
                        .HasColumnType("integer");

                    NpgsqlPropertyBuilderExtensions.UseIdentityByDefaultColumn(b.Property<int>("Id"));

                    b.Property<double>("Delta")
                        .HasColumnType("double precision");

                    b.Property<double>("Gamma")
                        .HasColumnType("double precision");

                    b.Property<double>("Rho")
                        .HasColumnType("double precision");

                    b.Property<double>("Theta")
                        .HasColumnType("double precision");

                    b.Property<double>("Unrealized_Pnl")
                        .HasColumnType("double precision");

                    b.Property<double>("Vega")
                        .HasColumnType("double precision");

                    b.HasKey("Id");

                    b.ToTable("OptionTradeEvaluation");
                });

            modelBuilder.Entity("Homework6.Trade", b =>
                {
                    b.Property<int>("Id")
                        .ValueGeneratedOnAdd()
                        .HasColumnType("integer");

                    NpgsqlPropertyBuilderExtensions.UseIdentityByDefaultColumn(b.Property<int>("Id"));

                    b.Property<int>("financialinstrumentid")
                        .HasColumnType("integer");

                    b.Property<double>("quantity")
                        .HasColumnType("double precision");

                    b.Property<double>("trade_price")
                        .HasColumnType("double precision");

                    b.HasKey("Id");

                    b.ToTable("Trade");
                });

            modelBuilder.Entity("Homework6.Unit", b =>
                {
                    b.Property<int>("Id")
                        .ValueGeneratedOnAdd()
                        .HasColumnType("integer");

                    NpgsqlPropertyBuilderExtensions.UseIdentityByDefaultColumn(b.Property<int>("Id"));

                    b.Property<double>("sizeUnit")
                        .HasColumnType("double precision");

                    b.Property<string>("typeUnit")
                        .IsRequired()
                        .HasColumnType("text");

                    b.HasKey("Id");

                    b.ToTable("Unit");
                });

            modelBuilder.Entity("Homework6.Option", b =>
                {
                    b.HasBaseType("Homework6.FinancialInstrument");

                    b.Property<DateTime>("expiration_date")
                        .HasColumnType("timestamp with time zone");

                    b.Property<int>("underlyingid")
                        .HasColumnType("integer");

                    b.HasIndex("underlyingid");

                    b.ToTable("Option");
                });

            modelBuilder.Entity("Homework6.Underlying", b =>
                {
                    b.HasBaseType("Homework6.FinancialInstrument");

                    b.Property<string>("Name")
                        .IsRequired()
                        .HasColumnType("text");

                    b.Property<string>("Symbol")
                        .IsRequired()
                        .HasColumnType("text");

                    b.ToTable("Underlying");
                });

            modelBuilder.Entity("Homework6.Asian", b =>
                {
                    b.HasBaseType("Homework6.Option");

                    b.Property<bool>("Is_Call")
                        .HasColumnType("boolean");

                    b.Property<double>("Strike")
                        .HasColumnType("double precision");

                    b.ToTable("Asian");
                });

            modelBuilder.Entity("Homework6.Barrier", b =>
                {
                    b.HasBaseType("Homework6.Option");

                    b.Property<double>("Barrier_Level")
                        .HasColumnType("double precision");

                    b.Property<bool>("Is_Call")
                        .HasColumnType("boolean");

                    b.Property<int>("Knock_Type")
                        .HasColumnType("integer");

                    b.Property<double>("Strike")
                        .HasColumnType("double precision");

                    b.ToTable("Barrier");
                });

            modelBuilder.Entity("Homework6.Digital", b =>
                {
                    b.HasBaseType("Homework6.Option");

                    b.Property<bool>("Is_Call")
                        .HasColumnType("boolean");

                    b.Property<double>("Payout")
                        .HasColumnType("double precision");

                    b.Property<double>("Strike")
                        .HasColumnType("double precision");

                    b.ToTable("Digital");
                });

            modelBuilder.Entity("Homework6.European", b =>
                {
                    b.HasBaseType("Homework6.Option");

                    b.Property<bool>("Is_Call")
                        .HasColumnType("boolean");

                    b.Property<double>("Strike")
                        .HasColumnType("double precision");

                    b.ToTable("European");
                });

            modelBuilder.Entity("Homework6.Lookback", b =>
                {
                    b.HasBaseType("Homework6.Option");

                    b.Property<bool>("Is_Call")
                        .HasColumnType("boolean");

                    b.Property<double>("Strike")
                        .HasColumnType("double precision");

                    b.ToTable("Lookback");
                });

            modelBuilder.Entity("Homework6.Range", b =>
                {
                    b.HasBaseType("Homework6.Option");

                    b.ToTable("Range");
                });

            modelBuilder.Entity("Homework6.FinancialInstrument", b =>
                {
                    b.HasOne("Homework6.Market", "MyMarket")
                        .WithMany()
                        .HasForeignKey("marketid")
                        .OnDelete(DeleteBehavior.Cascade)
                        .IsRequired();

                    b.Navigation("MyMarket");
                });

            modelBuilder.Entity("Homework6.Market", b =>
                {
                    b.HasOne("Homework6.Exchange", "MyExchange")
                        .WithMany()
                        .HasForeignKey("ExchangeId")
                        .OnDelete(DeleteBehavior.Cascade)
                        .IsRequired();

                    b.HasOne("Homework6.Unit", "MyUnit")
                        .WithMany()
                        .HasForeignKey("UnitId")
                        .OnDelete(DeleteBehavior.Cascade)
                        .IsRequired();

                    b.Navigation("MyExchange");

                    b.Navigation("MyUnit");
                });

            modelBuilder.Entity("Homework6.Option", b =>
                {
                    b.HasOne("Homework6.FinancialInstrument", null)
                        .WithOne()
                        .HasForeignKey("Homework6.Option", "Id")
                        .OnDelete(DeleteBehavior.Cascade)
                        .IsRequired();

                    b.HasOne("Homework6.Underlying", "underlying")
                        .WithMany()
                        .HasForeignKey("underlyingid")
                        .OnDelete(DeleteBehavior.Cascade)
                        .IsRequired();

                    b.Navigation("underlying");
                });

            modelBuilder.Entity("Homework6.Underlying", b =>
                {
                    b.HasOne("Homework6.FinancialInstrument", null)
                        .WithOne()
                        .HasForeignKey("Homework6.Underlying", "Id")
                        .OnDelete(DeleteBehavior.Cascade)
                        .IsRequired();
                });

            modelBuilder.Entity("Homework6.Asian", b =>
                {
                    b.HasOne("Homework6.Option", null)
                        .WithOne()
                        .HasForeignKey("Homework6.Asian", "Id")
                        .OnDelete(DeleteBehavior.Cascade)
                        .IsRequired();
                });

            modelBuilder.Entity("Homework6.Barrier", b =>
                {
                    b.HasOne("Homework6.Option", null)
                        .WithOne()
                        .HasForeignKey("Homework6.Barrier", "Id")
                        .OnDelete(DeleteBehavior.Cascade)
                        .IsRequired();
                });

            modelBuilder.Entity("Homework6.Digital", b =>
                {
                    b.HasOne("Homework6.Option", null)
                        .WithOne()
                        .HasForeignKey("Homework6.Digital", "Id")
                        .OnDelete(DeleteBehavior.Cascade)
                        .IsRequired();
                });

            modelBuilder.Entity("Homework6.European", b =>
                {
                    b.HasOne("Homework6.Option", null)
                        .WithOne()
                        .HasForeignKey("Homework6.European", "Id")
                        .OnDelete(DeleteBehavior.Cascade)
                        .IsRequired();
                });

            modelBuilder.Entity("Homework6.Lookback", b =>
                {
                    b.HasOne("Homework6.Option", null)
                        .WithOne()
                        .HasForeignKey("Homework6.Lookback", "Id")
                        .OnDelete(DeleteBehavior.Cascade)
                        .IsRequired();
                });

            modelBuilder.Entity("Homework6.Range", b =>
                {
                    b.HasOne("Homework6.Option", null)
                        .WithOne()
                        .HasForeignKey("Homework6.Range", "Id")
                        .OnDelete(DeleteBehavior.Cascade)
                        .IsRequired();
                });
#pragma warning restore 612, 618
        }
    }
}
